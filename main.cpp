#include <iostream>
#include <bits/stl_algo.h>
#include <sstream>

using namespace std;

typedef unsigned short elem_t;
const size_t elemSize = sizeof(elem_t);
const size_t B = 2;
const size_t B_BYTES = B * elemSize;
const size_t M_BYTES = 1 * 1024 * 1024 * 1024; //1GB
const size_t M = M_BYTES / elemSize;
const size_t m = M / B;


namespace util {
    double getTime() {
        struct timespec tp1;
        clock_gettime(CLOCK_REALTIME, &tp1);
        double theseSecs = tp1.tv_sec + tp1.tv_nsec / 1e9;
        return theseSecs * 1000;
    }
}

class BigFile {
public:
    static BigFile *getTempFile() {
        string name;
        stringstream strstream;
        strstream << "temp_" << tempfileCounter;
        strstream << ".bin";
        strstream >> name;
        tempfileCounter++;
        return new BigFile(name, true);
    }

    BigFile(string name, bool isTemp) : filename(name), isTemp(isTemp) {
        bool exists = false;
        if (FILE *test = fopen(filename.data(), "r")) {
            fclose(test);
            exists = true;
        }
        file = fopen64(filename.data(), exists ? "r+b" : "w+b");
    }

    virtual ~BigFile() {
        fclose(file);
        if (isTemp) {
            remove(filename.data());
        }
    }

    void writeB(void *ptr, size_t size, uint64_t pos) {
        double startTime = util::getTime();
        fseeko64(file, pos * B_BYTES, SEEK_SET);
        fwrite(ptr, 1, size * elemSize, file);
        ioTime += util::getTime() - startTime;
    }

    void startWrite(uint64_t pos) {
        buffer = (elem_t *) malloc(B_BYTES);
        position = pos;
        buffer_size = 0;
    }

    void write(elem_t e) {
        buffer[buffer_size++] = e;
        if (buffer_size == B) {
            writeB(buffer, buffer_size, position++);
            buffer_size = 0;
        }
    }

    void finishWrite() {
        if (buffer_size > 0) {
            writeB(buffer, buffer_size, position);
        }
        free(buffer);
    }

    void startRead() {
        buffer = (elem_t *) malloc(B_BYTES);
        position = 0;
        buffer_size = 0;
        current_index = 0;
    }

    elem_t read() {
        if (buffer_size == 0) {
            buffer_size = readB(buffer, position++);
            current_index = 0;
        }
        if (buffer_size == 0) throw "Error";
        buffer_size--;
        return buffer[current_index++];
    }

    void finishRead() {
        free(buffer);
    }

    size_t readB(void *ptr, uint64_t pos) {
        double startTime = util::getTime();
        fseeko64(file, pos * B_BYTES, SEEK_SET);
        size_t readCount = fread(ptr, elemSize, B, file);
        ioTime += util::getTime() - startTime;
        return readCount;
    }

    uint64_t readMore(void *ptr, uint64_t pos, size_t count) {
        uint64_t readElems = 0;
        for (size_t i = 0; i < count; i++) {
            size_t read_count = readB(ptr + readElems * elemSize, pos++);
            if (read_count == 0) {
                break;
            }
            readElems += read_count;
        }
        return readElems;
    }

    uint64_t get_N() {
        fseeko64(file, 0, SEEK_END);
        __off64_t size = ftello64(file);
        if (size < 0) throw "Error";
        return (uint64_t) (size / elemSize);
    }

    uint64_t get_n() {
        uint64_t N = get_N();
        uint64_t n = N / B;
        if (N % B) {
            n += 1;
        }
        return n;
    }

    void printFile(size_t maxItems) {
        elem_t *b = (elem_t *) malloc(B_BYTES);
        uint64_t blockCount = get_n();
        size_t counter = 0;
        printf("[");
        for (uint64_t i = 0; i < blockCount; i++) {
            size_t readCount = readB(b, i);
            for (size_t j = 0; j < readCount; ++j) {
                cout << " " << b[j] << " ";
                counter++;
                if (maxItems > 0 && counter == maxItems) {
                    goto BREAK;
                }
            }
        }
        BREAK:
        printf("]\n");
        free(b);
    }

    double getTotalIoTime() const {
        return ioTime;
    }

    string getFilename() const {
        return filename;
    }

private:
    FILE *file;
    string filename;
    elem_t *buffer;
    uint64_t position;
    size_t buffer_size;
    size_t current_index;
    double ioTime = 0;
    bool isTemp = false;
    static long tempfileCounter;
};

long BigFile::tempfileCounter = 0;


BigFile *buildU(BigFile *file, size_t a) {
    BigFile *U = BigFile::getTempFile();
    elem_t *read_buf = (elem_t *) malloc(M_BYTES);
    elem_t *write_buf = (elem_t *) malloc(B_BYTES);
    uint64_t n = file->get_n();
    for (uint64_t j = 0; j <= n; j += m) {
        size_t read_count = 0;
        for (size_t i = 0; i < m; i++) {
            size_t read = file->readB(read_buf + read_count, j + i); // * elemSize ??
            if (read == 0) {
                break;
            }
            read_count += read;
        }

        size_t write_buf_size = 0;
        uint64_t blocks_written = 0;
        for (size_t i = a - 1; i < read_count; i += a) {
            write_buf[write_buf_size++] = read_buf[i];
            if (write_buf_size == B) {
                U->writeB(write_buf, write_buf_size, blocks_written++);
                write_buf_size = 0;
            }
        }
        if (write_buf_size > 0) {
            U->writeB(write_buf, write_buf_size, blocks_written);
        }
    }
    free(read_buf);
    free(write_buf);
    return U;
}

BigFile **partition(BigFile *file, elem_t *pivots, size_t pivots_count) {
    size_t block_count = pivots_count + 1;
    BigFile **files = (BigFile **) malloc(block_count * sizeof(BigFile *));
    for (size_t j = 0; j < block_count; j++) {
        files[j] = BigFile::getTempFile();
        files[j]->startWrite(0);
    }
    elem_t *buffers = (elem_t *) malloc(B_BYTES * block_count);
    uint64_t N = file->get_N();
    file->startRead();
    for (uint64_t i = 0; i < N; i++) {
        elem_t e = file->read();
        for (size_t j = 0; j < block_count; j++) {
            if (j == pivots_count || e < pivots[j]) {
                files[j]->write(e); //TODO improve logic for e == pivot
                break;
            }
        }
    }
    file->finishRead();
    for (size_t j = 0; j < block_count; j++) {
        files[j]->finishWrite();
    }
    free(buffers);
    return files;
}

elem_t median_of_medians(BigFile *file) {
    uint64_t n = file->get_n();
    if (n < m) {
        elem_t *m_buf = (elem_t *) malloc(M_BYTES);
        uint64_t read_elems = file->readMore(m_buf, 0, (size_t) n);
        nth_element(m_buf, m_buf + read_elems / 2, m_buf + read_elems);
        elem_t result = m_buf[read_elems / 2];
        free(m_buf);
        return result;
    } else {
        elem_t *read_buf = (elem_t *) malloc(B_BYTES);
        BigFile *write_file = BigFile::getTempFile();
        write_file->startWrite(0);
        for (uint64_t i = 0; i < n; i++) {
            size_t readB = file->readB(read_buf, i);
            nth_element(read_buf, read_buf + readB / 2, read_buf + readB);
            elem_t median = read_buf[readB / 2];
            write_file->write(median);
        }
        write_file->finishWrite();
        elem_t result = median_of_medians(write_file);
        free(read_buf);
        delete write_file;
        return result;
    }
}

elem_t k_seq(BigFile *file, uint64_t k) {
    elem_t median = median_of_medians(file);
    elem_t pivot[] = {median};
    BigFile **files = partition(file, pivot, 1);
    BigFile *leftFile = files[0];
    BigFile *rightFile = files[1];
    free(files);
    uint64_t leftN = leftFile->get_N();
    uint64_t rightN = rightFile->get_N();
    elem_t kth;
    if (k <= leftN) {
        delete rightFile;
        kth = k_seq(leftFile, k);
    } else if (k > rightN) {
        delete leftFile;
        kth = k_seq(rightFile, k - leftN - 1);
    } else {
        delete leftFile;
        delete rightFile;
        kth = median;
    }
    return kth;
}

elem_t *getPivots(BigFile *file, size_t mu) {
    elem_t *pivots = (elem_t *) malloc(mu * elemSize);
    BigFile *U = buildU(file, sqrt(m) / 4);
    uint64_t N = U->get_N();
    uint64_t step = N / mu;
    for (size_t i = 0; i < mu; i++) {
        pivots[i] = k_seq(U, i * step);
    }
    delete U;
    return pivots;
}

void merge(BigFile **files, size_t count, BigFile *result) {
    result->startWrite(0);
    for (size_t i = 0; i < count; i++) {
        BigFile *file = files[i];
        uint64_t N = file->get_N();
        file->startRead();
        for (uint64_t j = 0; j < N; j++) {
            result->write(file->read());
        }
        file->finishRead();
    }
    result->finishWrite();
}

BigFile *distribution_sort(BigFile *file, BigFile *result) {
    size_t mu = sqrt(m);
    elem_t *pivots = getPivots(file, mu);
    BigFile **files = partition(file, pivots, mu);
    for (size_t i = 0; i <= mu; i++) {
        BigFile *temp = distribution_sort(files[i], NULL);
        delete files[i];
        files[i] = temp;
    }
    if (result == NULL) {
        result = BigFile::getTempFile();
    }
    merge(files, mu + 1, result);
    for (size_t i = 0; i <= mu; i++) {
        delete files[i];
    }
    free(files);
    free(pivots);
    return result;
}

void writeFile(const char* filename, uint64_t length) {
    BigFile* file = new BigFile(filename, false);
    file->startWrite(0);
    srand(100500);
    for (uint64_t i = 0; i < length; i++) {
        file->write((elem_t) rand());
    }
    file->finishWrite();
    delete file;
}

void test() {
    double startTime = util::getTime();
    BigFile file("int.bin", false);
    printf("block count:%lld\n", file.get_n());
    elem_t *b = (elem_t *) malloc(B_BYTES);
    size_t blockSize = file.readB(b, 2);
    printf("blocksize:%d, block=[", blockSize);
    for (size_t i = 0; i < blockSize; i++) {
        cout << " " << b[i] << " ";
    }
    printf("]\n");
    file.writeB(b, blockSize, 4);
    free(b);
    double allTime = util::getTime() - startTime;
    double ioTime = file.getTotalIoTime();
    double percentIO = ioTime / allTime * 100;
    printf("Total time:%f, IO: %f (%f%%)\n", allTime, ioTime, percentIO);
}


int main(int argc, char* argv[]) {
//    BigFile* file = new BigFile(argv[1], false);
//    BigFile* result = new BigFile(argv[2], false);
//    distribution_sort(file, result);
    writeFile("input.bin", 10);
    BigFile* file = new BigFile("input.bin", false);
    BigFile* result = new BigFile("result.bin", false);

    double startTime = util::getTime();
    distribution_sort(file, result);
    double allTime = util::getTime() - startTime;
    double ioTime = result->getTotalIoTime();
    double percentIO = ioTime / allTime * 100;
    file->printFile(0);
    result->printFile(0);
    printf("Total time:%f, IO: %f (%f%%)\n", allTime, ioTime, percentIO);
    return 0;
}
