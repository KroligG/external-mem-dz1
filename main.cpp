#include <iostream>
#include <bits/stl_algo.h>
#include <sstream>
#include <string.h>

using namespace std;

typedef long long elem_t;
const size_t elemSize = sizeof(elem_t);
const size_t B = 10000000;
const size_t B_BYTES = B * elemSize;
const size_t M_BYTES = 1 * 1024 * 1024 * 1536; //1.5GB
const size_t M = M_BYTES / elemSize;
const size_t m = M / B;
const size_t PARTITION_COUNT = 16;
const size_t P = 1000;
const char* TEMP_PATH;
//const size_t m = 144;
//const size_t M = m * B;
//const size_t M_BYTES = M * elemSize;


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
        strstream << TEMP_PATH << "temp_" << tempfileCounter;
        strstream << ".bin";
        strstream >> name;
        tempfileCounter++;
        return new BigFile(name, true);
    }

    BigFile(string name, bool isTemp) : filename(name), isTemp(isTemp) {
        bool exists = false;
        if (FILE *test = fopen64(filename.data(), "r")) {
            fclose(test);
            exists = true;
        }
        file = fopen64(filename.data(), exists ? "r+b" : "w+b");
        if (!file) {
            cout << strerror(errno) << endl;
            throw errno;
        }
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
        ioTime += (util::getTime() - startTime);
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
        ioTime += (util::getTime() - startTime);
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
        startRead();
        size_t readTo = (size_t) (maxItems ? maxItems : get_N());
        printf("[");
        for(size_t i = 0; i < readTo; i++) {
            cout << " " << read() << " ";
        }
        printf("]\n");
        finishRead();
    }

    void addAccumuatedIoTime(double t) {
        accumulatedIoTime += t;
    }

    double getTotalIoTime() const {
        return ioTime + accumulatedIoTime;
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
    double accumulatedIoTime = 0;
    bool isTemp = false;
    static long tempfileCounter;
};

long BigFile::tempfileCounter = 0;


BigFile *buildU(BigFile *file, size_t a) {
    BigFile *U = BigFile::getTempFile();
    uint64_t n = file->get_n();
    printf("Building U from file (block count %lld), taking every %d element\n", n, a);
    elem_t *m_buf = (elem_t *) malloc(M_BYTES);
    U->startWrite(0);
    for (uint64_t j = 0; j <= n; j += m) {
        uint64_t read_elems = file->readMore(m_buf, j, m);
        sort(m_buf, m_buf + read_elems);
        for (uint64_t i = a - 1; i < read_elems; i += a) {
            U->write(m_buf[i]);
        }
    }
    U->finishWrite();
    free(m_buf);
    return U;
}

BigFile **partition(BigFile *file, elem_t *pivots, size_t pivots_count, bool removePivots) {
    size_t block_count = pivots_count + 1;
    uint64_t N = file->get_N();
    printf("Partition file (elements %lld) to %d parts\n", N, block_count);
    BigFile **files = (BigFile **) malloc(block_count * sizeof(BigFile *));
    for (size_t j = 0; j < block_count; j++) {
        files[j] = BigFile::getTempFile();
        files[j]->startWrite(0);
    }
    file->startRead();
    int *pivotEqualitySwitches = (int *) malloc(pivots_count * sizeof(int));
    memset(pivotEqualitySwitches, 0 , pivots_count * sizeof(int));
    for (uint64_t i = 0; i < N; i++) {
        elem_t e = file->read();
        for (size_t j = 0; j < block_count; j++) {
            if (j == pivots_count || e < pivots[j]) {
                files[j]->write(e);
                break;
            }
            if(e == pivots[j]) {
                if(removePivots) {
                    break;
                }
                int swtch = pivotEqualitySwitches[j];
                files[j + swtch]->write(e);
                pivotEqualitySwitches[j] = swtch == 0 ? 1 : 0;
                break;
            }
        }
    }
    free(pivotEqualitySwitches);
    file->finishRead();
    for (size_t j = 0; j < block_count; j++) {
        files[j]->finishWrite();
    }
    return files;
}

elem_t median_of_medians(BigFile *file) {
    uint64_t n = file->get_n();
    printf("Find median of medians (file blocks %lld)\n", n);
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
        free(read_buf);
        elem_t result = median_of_medians(write_file);
        file->addAccumuatedIoTime(write_file->getTotalIoTime());
        delete write_file;
        return result;
    }
}

elem_t k_seq(BigFile *file, uint64_t k) {
    elem_t median = median_of_medians(file);
    elem_t pivot[] = {median};
    BigFile **files = partition(file, pivot, 1, true);
    BigFile *leftFile = files[0];
    BigFile *rightFile = files[1];
    free(files);
    uint64_t leftN = leftFile->get_N();
    uint64_t rightN = rightFile->get_N();
    printf("Find k-seq statistic, k=%lld, (left elems %lld) (right elems %lld)\n", k, leftN, rightN);
    elem_t kth;
    if (k <= leftN) {
        file->addAccumuatedIoTime(rightFile->getTotalIoTime());
        delete rightFile;
        rightFile = NULL;
        if(leftN < M) {
            elem_t *m_buf = (elem_t *) malloc(M_BYTES);
            uint64_t read_elems = leftFile->readMore(m_buf, 0, (size_t) leftFile->get_n());
            nth_element(m_buf, m_buf + k - 1, m_buf + read_elems);
            kth = m_buf[k - 1];
            free(m_buf);
        } else {
            kth = k_seq(leftFile, k);
        }
    } else if (k > rightN) {
        file->addAccumuatedIoTime(leftFile->getTotalIoTime());
        delete leftFile;
        leftFile = NULL;
        if(rightN < M) {
            uint64_t new_k = k - leftN - 1;
            elem_t *m_buf = (elem_t *) malloc(M_BYTES);
            uint64_t read_elems = rightFile->readMore(m_buf, 0, (size_t) rightFile->get_n());
            nth_element(m_buf, m_buf + new_k - 1, m_buf + read_elems);
            kth = m_buf[new_k - 1];
            free(m_buf);
        } else {
            kth = k_seq(rightFile, k - leftN - 1);
        }
    } else {
        kth = median;
    }
    if (leftFile != NULL) {
        file->addAccumuatedIoTime(leftFile->getTotalIoTime());
        delete leftFile;
    }
    if (rightFile != NULL) {
        file->addAccumuatedIoTime(rightFile->getTotalIoTime());
        delete rightFile;
    }
    return kth;
}

elem_t *getPivots(BigFile *file, size_t mu) {
    printf("Finding %d pivots\n", mu);
    elem_t *pivots = (elem_t *) malloc(mu * elemSize);
    BigFile *U = buildU(file, P);
    uint64_t N = U->get_N();
    uint64_t step = N / (mu + 1);
    for (size_t i = 1; i <= mu; i++) {
        pivots[i - 1] = k_seq(U, i * step);
    }
    file->addAccumuatedIoTime(U->getTotalIoTime());
    delete U;
    return pivots;
}

void merge(BigFile **files, size_t count, BigFile *result) {
    printf("Merging %d files\n", count);
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
    uint64_t n = file->get_n();
    printf("Distribution sort iteration, file blocks %lld\n", n);
    if (result == NULL) {
        result = BigFile::getTempFile();
    }
    if (n < m) {
        elem_t *m_buf = (elem_t *) malloc(M_BYTES);
        uint64_t read_elems = file->readMore(m_buf, 0, (size_t) n);
        sort(m_buf, m_buf + read_elems);
        result->startWrite(0);
        for (uint64_t i = 0; i < read_elems; i++) {
            result->write(m_buf[i]);
        }
        result->finishWrite();
        free(m_buf);
        return result;
    } else {
        size_t mu = PARTITION_COUNT - 1;
        elem_t *pivots = getPivots(file, mu);
        BigFile **files = partition(file, pivots, mu, false);
        free(pivots);
        for (size_t i = 0; i <= mu; i++) {
            BigFile *temp = distribution_sort(files[i], NULL);
            file->addAccumuatedIoTime(files[i]->getTotalIoTime());
            delete files[i];
            files[i] = temp;
        }
        merge(files, mu + 1, result);
        for (size_t i = 0; i <= mu; i++) {
            file->addAccumuatedIoTime(files[i]->getTotalIoTime());
            delete files[i];
        }
        free(files);
        return result;
    }
}

void writeFile(const char *filename, uint64_t length) {
    double startTime = util::getTime();
    BigFile *file = new BigFile(filename, false);
    file->startWrite(0);
    srand(100500);
    for (uint64_t i = 0; i < length; i++) {
//        file->write((elem_t) rand());
//        file->write((elem_t) (length - i));
        file->write((elem_t) rand() % length);
    }
    file->finishWrite();
    double allTime = util::getTime() - startTime;
    double ioTime = file->getTotalIoTime();
    double percentIO = ioTime / allTime * 100;
    double speed = ((double) (length * elemSize)) / (1024.0 * 1024.0) / (ioTime / 1000.0);
    printf("Total time:%f, IO: %f (%f%%), speed: %f Mb/s\n", allTime, ioTime, percentIO, speed);
    delete file;
}

void testFileReadSpeed(const char *filename) {
    double startTime = util::getTime();
    BigFile *file = new BigFile(filename, false);
    elem_t *m_buf = (elem_t *) malloc(M_BYTES);
    uint64_t read_elems = file->readMore(m_buf, 0, m);
    sort(m_buf, m_buf + read_elems);
    free(m_buf);
    double allTime = util::getTime() - startTime;
    double ioTime = file->getTotalIoTime();
    double percentIO = ioTime / allTime * 100;
    double speed = ((double) (read_elems * elemSize)) / (1024.0 * 1024.0) / (ioTime / 1000.0);
    printf("Total time:%f, IO: %f (%f%%), speed: %f Mb/s\n", allTime, ioTime, percentIO, speed);
    delete file;
}

void run_sort(const char* inputName, const char* resultName) {
    BigFile* file = new BigFile(inputName, false);
    BigFile* result = new BigFile(resultName, false);

    double startTime = util::getTime();

    distribution_sort(file, result);

    double allTime = util::getTime() - startTime;
    result->addAccumuatedIoTime(file->getTotalIoTime());
    double ioTime = result->getTotalIoTime();
    double percentIO = ioTime / allTime * 100;

    file->printFile(100);
    result->printFile(100);
    printf("Total time:%f, IO: %f (%f%%)\n", allTime, ioTime, percentIO);

    FILE *result_log = fopen("result_log.txt", "a");
    fprintf(result_log, "Total time: %f, IO: %f (%f%%), M: %d, B: %d, part: %d, P: %d\n", allTime, ioTime, percentIO, M, B, PARTITION_COUNT, P);
    fclose(result_log);
}

int main(int argc, char *argv[]) {
//    TEMP_PATH = "/media/lg/temp/";
    TEMP_PATH = argv[3];
//    writeFile("/media/lg/input2.bin", 1e8);
//    testFileReadSpeed("/media/lg/input.bin");
    run_sort(argv[1], argv[2]);
//    run_sort("/media/lg/input.bin", "/media/lg/result.bin");
//    run_sort("input.bin", "/media/lg/result2.bin");

    return 0;
}
