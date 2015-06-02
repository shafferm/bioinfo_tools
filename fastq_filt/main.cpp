#include <fstream>
#include <set>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

struct fastq_read {
    string id;
	CharString header;
	Dna5String seq;
	CharString qual;
} read1;

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
	if (argc!=4) {
		cout << "requires 3 arguements (Reads To_remove Output)" << endl;
		return -1;
	}

    SeqFileIn reader1(argv[1]);
    fstream in (argv[2]);
    SeqFileOut out1(argv[3]);
    
    // read in headers of files to toss
    std::set<string> headers;
    string line;
    while ( getline (in, line) ) {
        headers.insert(line);
    }
	
	while (!atEnd(reader1)) {
		readRecord(read1.header, read1.seq, read1.qual, reader1);
        read1.id = split(string(toCString(read1.header)), ' ')[0];
        if ( headers.end() == headers.find( read1.id ) )
            writeRecord(out1, read1.header, read1.seq, read1.qual);
    }    
}
