#ifndef RMD_H
#define RMD_H
#include <stdint.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <string.h>
#include <array>

using namespace std;

class WordBasedText;

static const int max_digit=40;

// Universal code
class UniversalCode {
public:
    UniversalCode(){buffer = (unsigned char*)calloc(50000000,1);buffer32=(uint32_t*)buffer;};
    virtual int encode(vector<uint64_t>){return 0;};
    virtual void code_init(vector<uint8_t>,int){};
    vector<uint32_t> code;
    uint32_t out[100];
    virtual void flush_to_byte(uint32_t){};
    int cur_byte=0,cur_bit=7,cur_value=0;
    virtual void reset() {cur_byte=0;cur_bit=7;curCodeSz=0;};
    virtual int process8bytes(int){return 0;};
    unsigned char* get_out(){return buffer;};
    virtual void dump(){};
    virtual int code_size(){return 0;};
    virtual uint32_t get_symbol(){return 0;};
    virtual void buidTableDecode(){};
    virtual ~UniversalCode(){delete[] buffer;};
    int serialize(uint8_t* external_buf){memcpy(external_buf,buffer,code_size());    return code_size();};
    void load(uint8_t* external_buf,int csize){memcpy(buffer,external_buf,csize);};
    virtual void codeIncrease(){};
    virtual uint8_t getDigitsNum(){return 0;};
    uint32_t* codewords;
    uint8_t* lengths;
    virtual void update_code(uint32_t*, uint32_t){};
    virtual void update_decode(uint32_t*, uint32_t){};
    virtual void print_code(){};
    virtual void initDecode(){};
/*private:*/
    unsigned char* buffer;
    uint32_t *buffer32,*decode;
    uint8_t *decodeL;
    uint32_t curCodeSz = 0;
    uint8_t maxClen;
protected:
    uint32_t sigma;
    bool smoothed_code=0;
};

class RMD: public UniversalCode
{
    public:
        RMD(vector<int>,int);
        virtual ~RMD();
        int encode(vector<uint64_t>);
        void code_output(int);
        void rmd_to_file(string);
        int code_size(){return cur_byte*4+1;};
    /*protected:

    private:*/
        int t=sizeof(unsigned int)*8;
        void flush_to_byte(uint32_t);
        void flush_to_byte32(uint32_t);
};

const int Tbyte_size=12,loop_n=64/Tbyte_size,Tsize=1<<Tbyte_size,Tshort_size=4,Tshort=1<<Tshort_size,Tmask=Tsize-1;

//===================== BCMix =============================
#define blockLen (10)
#define blockSz (1 << blockLen)
#define blockMask (blockSz - 1)

struct TableDecode
{
	uint32_t L[blockSz];
	uint8_t n[blockSz];
	uint8_t shift[blockSz];

	TableDecode* nxtTable;
};

// Binary mixed-digit code
class BCMix: public UniversalCode {
public:
    BCMix(vector<uint8_t>,int);
    virtual ~BCMix(){/*delete[] myCode; delete[] myCodeSz;*/};
    void reset() {pos=curCode=curSz=bitLen=bitStream=file8Pos=curCodeSz=0;};
    int code_size(){return (pos+1)*4;};
    uint32_t get_symbol();
    void code_init(vector<uint8_t>,int);
    void buidTableDecode();
    void codeIncrease();
    uint8_t getDigitsNum(){return DigitsNum[curCodeSz];};
    void flush_to_byte(uint32_t);

private:
    TableDecode ts[3];
    const static int max_digit=32;
    uint8_t blocksSz[max_digit];
    uint32_t pws[max_digit+10];
    vector<uint8_t> DigitsNum;
/*    uint32_t* myCode;
    uint32_t* myCodeSz;*/
    uint32_t pos = 0,file8Pos=0,bestBlocksNum;
    uint64_t curCode = 0;
    uint32_t curSz = 0;
    uint64_t bitStream = 0;
	uint32_t bitLen = 0;
	uint8_t prevCodeDigits = 0;

    void flush_to_byte32(uint32_t);
    void precalcCodesAndSizes(uint32_t);
};

class Shannon: public UniversalCode {
public:
    Shannon(uint32_t*,uint32_t,uint32_t);
    Shannon(uint32_t);
    Shannon(uint32_t,uint32_t);
    ~Shannon(){delete[] codewords; delete[] lengths; delete[] decodeL; delete[] decode;};
    void initCode();
    void print_code();
    void reset(){cur_bit=0;};
    int code_size(){return (cur_bit + 7) / 8;};
    void update_code(uint32_t*, uint32_t);
    void update_decode(uint32_t*, uint32_t);
//    uint32_t mask;
private:
};


#define CODES_PER_BITSTREAM 3
#define LOOKUP_BITS 12
#define LOOKUP_SIZE (1 << LOOKUP_BITS)
#define LOOKUP_MASK (LOOKUP_SIZE - 1)

class CHuffman: public UniversalCode {
public:
    /*void huffEncode(uint32_t* codes, uint32_t cntCodes, uint32_t dictSize, uint32_t* stream,
		uint32_t& streamLen, uint8_t* huffTable, uint32_t& huffTableLen);*/
	void precompute(uint8_t* huffTable, uint32_t huffTableLen, uint32_t dictSize);

	void huffDecode(uint32_t* codes, uint32_t cntCodes, uint32_t dictSize, uint32_t* stream,
		uint32_t streamLen);

    //Shannon(uint32_t*,uint32_t,uint32_t);
    CHuffman(uint32_t,int*);
    //Shannon(uint32_t,uint32_t);
    ~CHuffman(){delete[] codewords; delete[] lengths; delete[] decodeL; delete[] decode;};
    void initCode(uint32_t);
    void print_code(){};
    void reset(){cur_bit=0;};
    int code_size(){return (streamPos+1)*4;};
    void update_code(int*);
    void update_decode(uint32_t*, uint32_t);
    void initDecode();
    void flush_to_byte(uint32_t);
//    uint32_t mask;
protected:
    vector <array<int32_t, 2> > cntPar;
    vector <uint32_t> huffLens;
    void fillHuffLens(int32_t, uint32_t);
    uint32_t t_minHuffLen;
    uint32_t t_cntCodesPerLen[32];
    uint32_t t_firstCode[32], t_firstHuffCode[32];
    uint32_t t_firstCodeDiff[32]; // = t_firstCode[i] - t_firstHuffCode[i] (2 operation -> 1 operation)
    uint64_t t_limit[32];
    int8_t t_first[LOOKUP_SIZE];
private:
    uint8_t bufferShift;
    uint32_t streamPos;
    uint64_t buffer;
};


#endif
