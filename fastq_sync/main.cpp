#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <sstream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

struct fastq_read {
    string id;
	CharString header;
	Dna5String seq;
	CharString qual;
} read1, read2;

vector<string> &split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

vector<string> split(const string &s, char delim) {
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}

int main(int argc, char * argv[])
{
	if (argc!=6) {
		cout << "requires 5 arguements (R1,R2,R1S,R2S,SS)" << endl;
		return -1;
	}

    SeqFileIn reader1(argv[1]);
    SeqFileIn reader2(argv[2]);
    
    SeqFileOut out1(argv[3]);
	SeqFileOut out2(argv[4]);
	SeqFileOut out3(argv[5]);
    
    queue<string> previd1;
    map<string, fastq_read> prev1;
    queue<string> previd2;
    map<string, fastq_read> prev2;
	
	while (!atEnd(reader1) & !atEnd(reader2)) {
		readRecord(read1.header, read1.seq, read1.qual, reader1);
        read1.id = split(string(toCString(read1.header)), ' ')[0];
        readRecord(read2.header, read2.seq, read2.qual, reader2);
        read2.id = split(string(toCString(read2.header)), ' ')[0];
		cout << read1.id << endl;
		cout << read2.id << endl;
        if (read1.id == read2.id) {
			cout << "reads match" << endl;
            writeRecord(out1, read1.header, read1.seq, read1.qual);
            writeRecord(out2, read2.header, read2.seq, read2.qual);
        } else {
			cout << "reads dont match" << endl;
            //if read1.id in prev2 then print pair and all before are singles
            map<string,fastq_read>::iterator read2oldit = prev2.find(read1.id);
            if (read2oldit != prev2.end()) {
                //print read1 read2old pair
				fastq_read read2old = read2oldit->second;
                writeRecord(out1, read1.header, read1.seq, read1.qual);
                writeRecord(out2, read2old.header, read2old.seq, read2old.qual);
                while (previd2.front()!=read1.id) {
                    fastq_read read2old2 = prev2[previd2.front()];
                    writeRecord(out3, read2old2.header, read2old2.seq, read2old2.qual);
                    previd2.pop();
                    prev2.erase(read2old2.id);
                }
				previd2.pop();
				prev2.erase(read1.id);
            } else {
                previd1.push(read1.id);
                prev1.insert(map<string, fastq_read>::value_type(read1.id, read1));
            }
            map<string,fastq_read>::iterator read1oldit = prev1.find(read2.id);
            if (read1oldit != prev1.end()) {
				fastq_read read1old = read1oldit->second;
                writeRecord(out1, read1old.header, read1old.seq, read1old.qual);
                writeRecord(out2, read2.header, read2.seq, read2.qual);
                while (previd1.front()!=read2.id) {
                    fastq_read read1old1 = prev1[previd1.front()];
                    writeRecord(out3, read1old1.header, read1old1.seq, read1old1.qual);
                    previd1.pop();
                    prev1.erase(read1old1.id);
                }
				previd1.pop();
				prev1.erase(read2.id);
            } else {
                previd2.push(read2.id);
                prev2.insert(map<string, fastq_read>::value_type(read2.id, read2));
            }
        }
    }
    //WRITE IN END CONDITIONS FOR BOTH FILES
    if (atEnd(reader1)) {
		cout << "extra reads in 2" << endl;
        while (!atEnd(reader2)) {
	        readRecord(read2.header, read2.seq, read2.qual, reader2);
	        read2.id = split(string(toCString(read2.header)), ' ')[0];
			map<string,fastq_read>::iterator read1oldit = prev1.find(read2.id);
            if (read1oldit != prev1.end()) {
				fastq_read read1old = read1oldit->second;
                writeRecord(out1, read1old.header, read1old.seq, read1old.qual);
                writeRecord(out2, read2.header, read2.seq, read2.qual);
                while (previd1.front()!=read2.id) {
                    fastq_read read1old1 = prev1[previd1.front()];
                    writeRecord(out3, read1old1.header, read1old1.seq, read1old1.qual);
                    previd1.pop();
                    prev1.erase(read1old1.id);
                }
				previd1.pop();
				prev1.erase(read2.id);
            } else {
                previd2.push(read2.id);
                prev2.insert(map<string, fastq_read>::value_type(read2.id, read2));
            }
        }
    } else {
		cout << "extra reads in 1" << endl;
    	while (!atEnd(reader1)) {
			readRecord(read1.header, read1.seq, read1.qual, reader1);
	        read1.id = split(string(toCString(read1.header)), ' ')[0];
            map<string,fastq_read>::iterator read2oldit = prev2.find(read1.id);
            if (read2oldit != prev2.end()) {
                //print read1 read2old pair
				fastq_read read2old = read2oldit->second;
                writeRecord(out1, read1.header, read1.seq, read1.qual);
                writeRecord(out2, read2old.header, read2old.seq, read2old.qual);
                while (previd2.front()!=read1.id) {
                    fastq_read read2old2 = prev2[previd2.front()];
                    writeRecord(out3, read2old2.header, read2old2.seq, read2old2.qual);
                    previd2.pop();
                    prev2.erase(read2old2.id);
                }
				previd2.pop();
				prev2.erase(read1.id);
            } else {
                previd1.push(read1.id);
                prev1.insert(map<string, fastq_read>::value_type(read1.id, read1));
            }
    	}
    }
	//write remaining entries in both queues to singles file
	while (!previd1.empty()) {
		fastq_read single = prev1[previd1.front()];
		writeRecord(out3, single.header, single.seq, single.qual);
		previd1.pop();
		prev1.erase(single.id);
	}
	while (!previd2.empty()) {
		fastq_read single = prev2[previd2.front()];
		writeRecord(out3, single.header, single.seq, single.qual);
		previd2.pop();
		prev2.erase(single.id);
	}
}
