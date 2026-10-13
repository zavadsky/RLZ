#ifndef WORDBASEDTEXT_H
#define WORDBASEDTEXT_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <math.h>

using namespace std;

class WordBasedText
{
    public:
        int glob=0;
        WordBasedText(string,char);
        virtual ~WordBasedText();
        int word_frequences();
        map<string,int> word_freq; // <word,frequency> map
        void text_rewind(){buffer.clear(std::stringstream::goodbit); buffer.seekg(0); pos=0;};
        bool eof(){return buffer.eof();};
        string get_word();
        char get_char(){char c; buffer>>noskipws>>c; /*cout<<txt[pos]<<" "<<int(txt[pos])<<"|";*/ return c;/*return txt[pos++];*/};
        int getMaxSymb(){return diff_words;};
        map<string,int> word_symbol; // map <word; number in the list, ordered by frequency>
        vector<uint64_t> Frequencies; // vector of descending frequencies of all words
        map<int,int> freq_freq; // <frequency, number of words of this frequency> map
        int Nwords,Nchar; //total number of words and characters in the text
        void output_stream(string);
        vector<int> numbers;
        void prepare_adaptive();
        stringstream buffer;
        void write_to_file();

    protected:
        multimap<int,string,greater<int>> freq_word; // <frequency,word> multimap
        vector<int> DiffFreq; // vector of distinct frequencies
        int NFreq;   // number of distinct frequencies
        int diff_words=0;
        char alpha_num; // # - punctuation; 0 - non-alphanumeric; other - alphanumeric;

    private:
        unsigned char* txt;
        double entropy;
        string reserve_str;
        uint32_t pos=0;
};



#endif // WORDBASEDTEXT_H
