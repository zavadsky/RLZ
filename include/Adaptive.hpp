#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <map>
#include <set>
#include <deque>
#include <list>
#include <ostream>
#include <unordered_map>
#include <vector>
#include "UniversalCode.h"
#include "WordBasedText.h"
#include "sparse_table.hpp"

using namespace std;

extern double m_time;
extern std::chrono::high_resolution_clock::time_point start,stop;

#define measure_start()  \
start = std::chrono::high_resolution_clock::now();

#define measure_end()     stop = std::chrono::high_resolution_clock::now(); \
    m_time = (double)chrono::duration_cast<std::chrono::nanoseconds>( stop - start ).count()/1000000;

typedef unsigned char uchar;
typedef unsigned long long ull;

class BackwardEncoder {
public:
    BackwardEncoder(WordBasedText*);
    BackwardEncoder(int);
    BackwardEncoder(WordBasedText*,vector<uint32_t>);
    vector<uint32_t> get_dic(){return dic;};
    void checkDecode();
    void checkDecodeChar(WordBasedText*);
    virtual void serial(string){};
    void serializeChar(uint8_t*);
    void loadChar(uint8_t*);
protected:
    uint32_t n; // text length in characters
    unsigned char* outChar = nullptr;
    uint32_t* codeStream = nullptr;
    WordBasedText* text;
    unsigned char* inputText = nullptr;
    unsigned NYT;
    vector<uint32_t> out;
    vector<uint32_t> dic;
    ~BackwardEncoder(){
        delete[] inputText;
        };
};

typedef map<string,uint32_t> MSI;

struct Arc {
    string label;
    bool extra;
};

typedef struct Arc arc;


class TextParser : public BackwardEncoder {
static const int krm = 10299827;
//static const int krm = 27;
//static const int krm = 9000013;
private:
    static const int depth = 14/*, checkL = 4, Branches = 8, step = 4*/;
    int eulertigs_k = 4;
    // Counts and fingerprints are separate so frequencies retain all 32 bits.
    std::unique_ptr<uint32_t[][krm]> KRF;
    std::unique_ptr<uint16_t[][krm]> KRF_fingerprints;
    std::unique_ptr<uint64_t[]> KRF_lang_bits;
    std::vector<uint64_t> overflow_keys[depth];
    std::vector<uint32_t> overflow_counts[depth];
    std::vector<uint8_t> overflow_in_dictionary[depth];
    static const int /*L = 200,*/ char_L=7, mask = (1<<char_L)-1;
    static const int kra = 31;
    static const int const_len = 65300;
    static const uint32_t undef = 0xffffffff;
    int len = const_len,extra_len = 0;
    //const static int max_digit=40;
    int *f = nullptr,*f1 = nullptr;
    int L,segn,max_len;
    bool blocked_dictionary_mode = false;
    size_t common_dictionary_size = 0;
    size_t regular_block_part_size = 0;
    int dictionary_block_count = 1;
    int text_block_size = 0;
    int current_text_end = 0;
    double blocked_sa_time_ms = 0.0;
    double blocked_lcp_time_ms = 0.0;
    double blocked_rmq_time_ms = 0.0;
    double blocked_parse_time_ms = 0.0;
    double dictionary_stat_time_ms = 0.0;
    double dictionary_growth_time_ms = 0.0;
    double liao_index_time_ms = 0.0;
    double liao_parse_time_ms = 0.0;
    int err_num;
    int prev_err = -1;
    double minl;
    int Nsymb;
    MSI stat[depth];
    map<ull,int> stat1[8];
    map<ull,map<ull,int>> stat8[depth-8];
    map<ull,uchar> st8;
    //int stat1[depth][256];
    map<string,int> imbalance;
    string lang = " the that an";
    //string lang = " th";
    // Dynamically sized dictionary. Up to 256 missing byte values can be
    // appended, and depth bytes of zero padding protect fixed-width reads.
    uchar* langstr;
    size_t lang_capacity;
    //uchar lang1[len+255];
    set<unsigned char> symbols;
    vector<int> err;
    int iter(int);
    int iter_ii(int);
    unsigned char *result,*output;
    vector<pair<unsigned char,double>> err_masks;
    vector<int> improve;
    int iter_log(int);
    multimap<string,arc> deBrujn;
    list<arc> cycle;
    uint32_t *Lengths;
    uint32_t *Offsets;

    int bestFullBitsSz = 10000000,bestBlocksNum;
    void BCS(uint8_t[], uint32_t, uint64_t, uint64_t, uint64_t, uint64_t);
    int iter_sa(int);
    int iter_sa1(int);
    int iter_sa2(int);
    int iter_sa3(int);
    vector<uint64_t> PSum;
    uint64_t  getSum(uint32_t, uint32_t);
    vector<uint8_t> bestDigitsSz = {9,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    uint16_t KRDegrees[depth];
    uint16_t KRV[depth];
    int KarpRabin(int,int,int,int);
    int lcp(uchar*,uchar*);
    int lcp1(uchar*,uchar*);
    int bounded_lcp(int lang_pos, int text_pos) const;
    void build_hash_stat_range(uint32_t begin, uint32_t end);
    void ensure_hash_storage();
    string grow_hash_dictionary(const string& initial, size_t target_size);
    void build_unresolved_hash_stat_range(uint32_t begin, uint32_t end);
    string grow_unresolved_hash_dictionary(const string& initial,
                                             size_t target_size);
    void process_text_lz_blocked(int repetitions, bool verbose);
    void set_symbol_positions(const uchar* dictionary, size_t size);
    void KarpRabinInit();
    set<string> LangStrings[depth];
    ull valj(ull,ull,uint8_t) const;
    ull valj_full(ull,ull,uint8_t) const;
    void calculate_suffix_hashes(const uchar* window, ull hashes[depth]) const;
    static ull mix_hash(ull value);
    void hash_cells(ull hash, size_t cells[4]) const;
    int find_hash_cell(int index, ull hash) const;
    uint32_t hash_count(int index, ull hash) const;
    void increment_hash_count(int index, ull hash);
    bool hash_is_in_dictionary(int index, ull hash) const;
    void mark_hash_in_dictionary(int index, ull hash);
    uint16_t substring_fingerprint(const uchar* data, int length) const;
    uint32_t substring_primary_bucket(const uchar* data, int length) const;
    uint32_t chained_count(int index, uint32_t bucket,
                           uint16_t fingerprint) const;
    bool chained_in_dictionary(int index, uint32_t bucket,
                               uint16_t fingerprint) const;
    void mark_chained_in_dictionary(int index, uint32_t bucket,
                                    uint16_t fingerprint);
    ull vall(ull,ull,uint8_t);
    int symbol_pos[256];
    int SA_pos[257];
    uint64_t SA_pos2[65537];
    int64_t *SA = nullptr;
    int64_t *LCP = nullptr;
    SparseTable<int64_t*> *RMQ = nullptr;
    bool mode = 1; // mode=1 => stat is built over whole text
public:
    void report_hash_collisions(std::ostream& output) const;
    uint32_t* lengths_data(){return Lengths;};
    uint32_t* offsets_data(){return Offsets;};
    const uint8_t* dictionary_data() const { return langstr; }
    const uint8_t* text_data() const { return inputText; }
    int  segments_number(){return segn;};
    int dictionary_size() const { return len; }
    bool has_block_dictionary() const { return blocked_dictionary_mode; }
    size_t common_dictionary_bytes() const { return common_dictionary_size; }
    size_t block_part_bytes() const { return regular_block_part_size; }
    int dictionary_blocks() const { return dictionary_block_count; }
    size_t logical_dictionary_bytes() const {
        if (!blocked_dictionary_mode || dictionary_block_count <= 1)
            return static_cast<size_t>(len);
        const size_t last_part = static_cast<size_t>(len) -
            common_dictionary_size - regular_block_part_size *
                static_cast<size_t>(dictionary_block_count - 1);
        return common_dictionary_size +
            std::max(regular_block_part_size, last_part);
    }
    double blocked_sa_time() const { return blocked_sa_time_ms; }
    double blocked_lcp_time() const { return blocked_lcp_time_ms; }
    double blocked_rmq_time() const { return blocked_rmq_time_ms; }
    double blocked_parse_time() const { return blocked_parse_time_ms; }
    double dictionary_stat_time() const { return dictionary_stat_time_ms; }
    double dictionary_growth_time() const {
        return dictionary_growth_time_ms;
    }
    double liao_index_time() const { return liao_index_time_ms; }
    double liao_parse_time() const { return liao_parse_time_ms; }
    void init_SA();
    void init_LCP();
    void init_RMQ();
    void init_SA_old();
    void Eulerise();
    void Euler();
    void process_text_log();
    void process_text2_log();
    void encode();
    void build_stat();
    void build_stat1();
    void make_lang();
    void make_lang1();
    void make_lang2();
    void build_statistical_dictionary(size_t dictionary_size);
    void build_exact_dictionary(size_t dictionary_size);
    void build_unresolved_hash_dictionary(size_t dictionary_size);
    void build_unresolved_hash_dictionary(size_t dictionary_size,
                                          size_t hash_modulus);
    void build_suffix_array_dictionary(size_t dictionary_size,
                                       uint32_t max_suffix_length = 24,
                                       int length_correction = -1,
                                       const std::string& seed = "the ");
    void build_eulertigs_dictionary(int order = 4);
    void set_dictionary(const std::vector<uint8_t>& dictionary);
    void build_block_dictionary(size_t common_size,
                                size_t logical_dictionary_size,
                                int block_count);
    void process_text();
    void process_text2();
    void process_text_lz(int, bool verbose = true);
    int process_text_lz_arrays_only(int repetitions = 1);
    void process_text_liao(int repetitions = 1, bool verbose = true);
    void encode_values();
    void completeLang();
    void serialize_values();
    void serializeLang(string);
    void deserializeLang(string);
    void print_map();
    void calc_imbalance();
    TextParser(WordBasedText*);
    // Construct a deterministic dictionary from k-mers whose starting
    // positions are uniformly distributed over the input text.
    TextParser(WordBasedText*, unsigned kmer_length,
               size_t target_dictionary_size);
    // Use a text prefix followed by every text byte absent from that prefix.
    // The final dictionary has exactly target_dictionary_size bytes.
    TextParser(WordBasedText*, size_t target_dictionary_size);
    TextParser(WordBasedText*,string);
    ~TextParser();
    bool dictionary_covers_text_alphabet() const;
    void precalcBestBlocks();
    bool check_substr(string);
    void Spell();
    int avg;
};

// Our method
class AdaptiveN : public BackwardEncoder {
protected:
    vector<int> m;
    vector<int> rev;
    vector<int> ind;
    vector<int> f,f1;
    UniversalCode* C;
    vector<uint8_t> bestDigitsSz;
    int j;
    const static int max_var_digit=12,init_j=26;
    void update1(int);
    void update_code(uint32_t);
    void update_decode(uint32_t);
    void update_code_char(uint32_t);
    void update_decode_char(uint32_t);
    void precalcBestBlocks();
    void BCS(uint8_t[], uint32_t, uint64_t, uint64_t, uint64_t, uint64_t);
    const static int max_digit=40;
    int bestFullBitsSz = 10000000,bestBlocksNum;
    static const int sigma = 128;
public:
    int mc=0;
    uint32_t* PSum;
    uint64_t  getSum(uint32_t, uint32_t);
    double rate=5;
    AdaptiveN(WordBasedText*,UniversalCode*);
    AdaptiveN(UniversalCode*,int);
    AdaptiveN(WordBasedText*,UniversalCode*,vector<uint32_t>);
    double encode();
    double encode1();
    double encode2();
    double encode2_char();
    void print_code();
    void update();
    int decode();
    int decode1();
    int decode2();
    ~AdaptiveN(){delete[] outChar;};
};

class AdaptiveNDecode: public AdaptiveN {
private:
    int n;
public:
    int decode2_char();
    AdaptiveNDecode(int,UniversalCode*);
    ~AdaptiveNDecode(){delete[] outChar;};
};

typedef struct {
    unsigned up,      // next node up the tree
        down,         // pair of down nodes
        symbol,       // node symbol value
        weight;       // node weight
} HTable;

typedef struct {
    unsigned esc,     // the current tree height
        root,         // the root of the tree
        size,         // the alphabet size
        *map;         // mapping for symbols to nodes
    HTable table[1];  // the coding table starts here
} HCoder;

// Classic Vitter algorithm
class Vitter: public BackwardEncoder {
protected:
    FILE *In, *Out;
    unsigned char ArcBit = 0;
    uint32_t ArcChar = 0;
    HCoder *huff;
    WordBasedText* text;
    int len=0;
    const int size = (1<<24);
    void update(string);
    HCoder *huff_init (unsigned size, unsigned root);
    void arc_put1 (unsigned bit);
    unsigned huff_split (HCoder *, unsigned);
    unsigned huff_leader (HCoder *, unsigned);
    unsigned huff_slide (HCoder *, unsigned);
    void huff_increment (HCoder *, unsigned);
    void huff_scale (HCoder *, unsigned);
    void huff_sendid (HCoder *, unsigned);
    void huff_encode (HCoder *, unsigned);
public:
    Vitter(WordBasedText*);
    Vitter(int);
    ~Vitter(){delete[] huff;};
    double encode();
    int encode_char();
};

class VitterDecode: public Vitter {
public:
    VitterDecode(int);
    ~VitterDecode(){delete[] outChar;};;
    unsigned huff_decode (HCoder*);
    void decode_char();
private:
    unsigned arc_get1();
    unsigned huff_readid (HCoder*);
    uint32_t outPos,inPos;
};

// Modified Vitter algorithm
class VitterM : public BackwardEncoder {
private:
    FILE *In, *Out;
    unsigned char ArcBit = 0;
    int ArcChar = 0;
    HCoder *huff;
    int len=0,mc=0;
    const int size = (1<<24);
    void update(string);
    HCoder *huff_init (unsigned size, unsigned root);
    void arc_put1 (unsigned bit);
    unsigned arc_get1 ();
    unsigned huff_split (HCoder *, unsigned);
    unsigned huff_leader (HCoder *, unsigned);
    unsigned huff_slide (HCoder *, unsigned);
    void huff_increment (HCoder *, unsigned);
    void huff_sendid (HCoder *, unsigned);
    void huff_encode (HCoder *, unsigned);
    unsigned huff_decode (HCoder *);
    unsigned huff_readid (HCoder*);
    unsigned char *buffer;
public:
    VitterM(WordBasedText*);
    double encode();
    unsigned decode ();
    int serialize(uint8_t*);
    void load(uint8_t*,int);
    ~VitterM(){delete huff;};
};


class Gagie : public BackwardEncoder {
private:
protected:
    char code_type; // 0 - Huffman, 1 - Huffman smoothed, 2 - Shannon, 3 - Shannon smoothed
    uint32_t *freqs;
    UniversalCode* C;
    uint8_t bufferShift;
    uint32_t streamPos;
    uint64_t buffer;
    void BlockEncode(int);
    static const int sigma = 128;
    uint32_t inPos;
public:
    Gagie(WordBasedText*,char);
    Gagie(int text_len):BackwardEncoder(text_len){};
    int encode_char();
    int encode_char1();
    void encode_word();
    double entropy();
    ~Gagie();
};

class GagieDecode : public Gagie {
private:
    uint32_t outPos;
    void BlockDecode(int);
public:
    GagieDecode(int,int);
    void decode_char();
    void decode_char1();
    double entropy();
    ~GagieDecode(){delete[] outChar;};
};

void print_vector(vector<uint32_t>,string);
void print_vector(uint32_t*,uint32_t);
