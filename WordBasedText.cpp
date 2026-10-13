#include "include/WordBasedText.h"

WordBasedText::WordBasedText(string fname,char t)
{
    FILE *f;
    f = fopen(fname.c_str(),"rb");
    std::ifstream file(fname,ios::binary);
    buffer << file.rdbuf();
    file.seekg(0, std::ios::end);
    Nchar = file.tellg();
    alpha_num=t;
    text_rewind();
}


bool isalnum0(char c,char prev) {
    return isalnum(c) || c=='-' || (int(c)==39 && isalnum(prev));
}

// Reading a word from alpha-numeric or punctuation stream
string WordBasedText::get_word() {
string sa="",sn="";
char a,c=0,prev;
    if(buffer.eof())
        return "";
    if(alpha_num==2) {
        buffer>>sa;
        return sa;
    }
    if(alpha_num=='#') {
        buffer>>noskipws>>c;
        if(c) {
        sn=c;
        while(c!='#' && !buffer.eof()) {
            buffer>>noskipws>>c;
            sn+=c;
        }
        }
        return sn;
    }
    do {
        prev=c;
        buffer>>noskipws>>c;
        sn+=c==10?'#':c;
    } while(!isalnum0(c,prev) && ! buffer.eof());
    if(buffer.eof())
        return sn;
    buffer.unget();
    sn=sn.substr(0,sn.length()-1);
    do {
        prev=c;
        buffer>>noskipws>>c;
        sa+=c;
    } while(isalnum0(c,prev) && ! buffer.eof());
    if(buffer.eof())
        return sa;
    buffer.unget();
    if(alpha_num) {
        sa=sa.substr(0,sa.length()-1);
        return sa;
    }
    if(sn!="")
        return sn;
    do {
        prev=c;
        buffer>>noskipws>>c;
        sn+=c==10?'#':c;
    } while(!isalnum0(c,prev) && ! buffer.eof());
    if(buffer.eof())
        return sn;
    buffer.unget();
    sn=sn.substr(0,sn.length()-1);
    return sn;
}

//Create different maps and dictionaries
int WordBasedText::word_frequences() {
double pi,sum1=0,sum2=0;
string word;
int i=0,k=0;
    Nwords=0;
    // Create the map <word,frequency> - word_freq
	while ( ! buffer.eof() ) {
		word=get_word();
		if(word_freq.find(word)!=word_freq.end()) {
			word_freq[word]++;
        } else {
			word_freq.insert(make_pair(word,1));
		}

		Nwords++;
    }
	cout<<endl<<"Input file processed. <word,frequency> map buit. Words in the text: "<<Nwords<<endl;
	// Create the multimap <frequency,word> consisting all different words - freq_word
	for(auto it=word_freq.begin();it!=word_freq.end();it++,diff_words++) {
		freq_word.insert(make_pair(it->second,it->first));
		// Calculate Shannon entropy
		pi=(double)it->second;
		entropy-=pi*(long double)log2((long double)pi/Nwords);
	}
    cout<<"<frequency,word> map built "<<sum1<<" "<<sum2<<endl;
	// Create 1) map <word,symbol> which maps words of text to integers according to descending order of their frequencies;
	// 2) vector of Frequencies of all unique words; 3) vector of different Frequencies DiffFreq; 4) map freq_freq<frequency, number of words having this frequency>
int j=-1;
int frq=0;
NFreq=0;
multimap<int,string> :: iterator it1;
	for(it1=freq_word.begin(),i=0;it1!=freq_word.end();i++,it1++) {
        if(it1->first!=frq) {
            if(frq)
                freq_freq.insert(make_pair(frq,i-j));
            frq=it1->first;
            j=i;
            DiffFreq.push_back(frq);
            NFreq++;
        }
		word_symbol.insert(make_pair(it1->second,i)); //insert first element from freq_word map
		Frequencies.push_back(it1->first);
	}
    freq_freq.insert(make_pair(frq,i-j));
	cout<<"<word,symbol> map built. Different words in text: "<<diff_words<<". Entropy H0="<<(int)entropy/8<<" bytes."<<endl;
	return Nwords;
}

// Convert words to numbers accordingly to their leftmost occurrences
void WordBasedText::prepare_adaptive() {
map <string,int> occ;
int mc=0;
    text_rewind();
	while(!eof()) {
        string s=get_word();
        if(occ.find(s)==occ.end()) {
            occ.insert(make_pair(s,mc));
            numbers.push_back(mc++);
        } else
            numbers.push_back(occ[s]);
	}
}

void WordBasedText::write_to_file() {
ofstream myfile;
  myfile.open ("words.txt");
  text_rewind();
  while(!eof()) {
        myfile<<get_word()<<" ";
	}
   myfile.close();
}

WordBasedText::~WordBasedText()
{
    //dtor
}

