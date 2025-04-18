#pragma once
#include <bits/stdc++.h>
using namespace std;
struct Tree
{
    int frequency;
    unsigned char character;
    Tree *left = NULL;
    Tree *right = NULL;
};

class TreeComparator
{
  public:
    bool operator()(Tree *a, Tree *b)
    {
        return a->frequency > b->frequency;
    }
};

Tree *buildHuffmanTree(vector<int> freqTable)
{

    priority_queue<Tree *, vector<Tree *>, TreeComparator> huffqueue;
    for (int i = 0; i < freqTable.size(); i++)
    {
        Tree *node = new Tree();
        node->frequency = freqTable[i];
        node->character = (unsigned char) i;
        
        huffqueue.push(node);
    }

    

    while (huffqueue.size() > 1)
    {
        Tree *a, *b;
        a = huffqueue.top();
        huffqueue.pop();
        
        b = huffqueue.top();
        huffqueue.pop();
        Tree *c = new Tree();
        c->frequency = a->frequency + b->frequency;
        c->left = a;
        c->right = b;
        huffqueue.push(c);
    }

    return huffqueue.top();
}

string toBinary(unsigned  char a)
{
    string output  = "";
    while(a!=0)
    {
        string bit = a%2==0?"0":"1";
        output+=bit;
        a/=2;
    }

    if(output.size()<8)
    {
        int deficit = 8 - output.size();
        for(int i=0; i<deficit; i++)
        {
            output+="0";
        }
    }

    reverse(output.begin(), output.end());
    return output;
    
}

void traverseHuffmanTree(Tree *root, string prev, string toAppend, map<unsigned char, string> &codemap)
{

    prev+=toAppend;
    
    if (root->right == NULL && root->left == NULL)
    {
        // cout<<root->character<<" "<<prev<<endl;   
        codemap[root->character] = prev;
    }
    if (root->right != NULL)
    {
        traverseHuffmanTree(root->right, prev, "1", codemap);
    }

    if (root->left != NULL)
    {
        traverseHuffmanTree(root->left, prev, "0", codemap);
    }
}

unsigned char *readFileIntoBuffer(char *path, int &sz)
{
    FILE *fp = fopen(path, "rb");
    sz = 0;
    fseek(fp, 0, SEEK_END);
    sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    unsigned char *buffer = (unsigned char *)malloc(sz);
    fread(buffer, 1, sz, fp);
	fclose(fp);
    return buffer;
}

void writeFileFromBuffer(char *path, unsigned char *buffer, int sz, int flag)
{
    FILE *fp;
    if(flag==0)
    {
        fp = fopen(path, "wb");
    }

    else{
        fp = fopen(path, "ab");
    }
   
    fwrite(buffer, 1, sz, fp);

    fclose(fp);
}


vector<pair<unsigned char, int> > convertToVector(map<unsigned char, int> codes)
{
    vector<pair<unsigned char, int> > codesV;

    for (map<unsigned char, int>::iterator i = codes.begin(); i != codes.end(); i++)
    {
        codesV.push_back(make_pair(i->first, i->second));
    }

    return codesV;
}

string getHuffmanBitstring(unsigned char *buffer, map<unsigned char, string> codes, int sz, int& paddedBits)
{
    string outputBuffer="";
    for(int i=0; i<sz; i++)
    {
        outputBuffer.append(codes[buffer[i]]);
    }

    if(outputBuffer.size()%8!=0)
    {
        int deficit = 8*((outputBuffer.size()/8)+1)-outputBuffer.size();
        paddedBits = deficit;
        for(int i=0; i<deficit; i++)
        {
            outputBuffer+="0";
        }
    }

    //Instead of adding zeroes to the end, add zero to the begining of the last byte

    return outputBuffer;
    
}

unsigned char* getBufferFromString(string bitstring, vector<unsigned char>&outputBuffer, int& sz)
{
    int interval = 0;
    unsigned char bit = 0;

    for(int i=0; i<sz; i++)
    {
         bit = (bit<<1)|(bitstring[i]-'0');
         
        interval++;
        if(interval==8)
        {
            interval = 0;
            outputBuffer.push_back(bit);
            bit = 0;
        
        }   
    }
    sz = outputBuffer.size();
    return outputBuffer.data();
}

string getStringFromBuffer(unsigned char* buffer, int sz)
{
    string bitstring = ""; 
    for(int i=0; i<sz; i++)
    {
        bitstring+=toBinary(buffer[i]);
    }

    return bitstring;
}

unsigned char* getDecodedBuffer(string bitstring, vector<unsigned char>&buffer, map<unsigned char, string> codes, int &sz, int paddedBits)
{
    string bit = "";
    map<string, unsigned char> reversecodes;
    
    for(map<unsigned char, string>::iterator i = codes.begin(); i!=codes.end(); i++)
    {
        reversecodes[i->second] = i->first;
    }

    for(int i=0; i<bitstring.size()-paddedBits; i++)
    {
        bit+=string(1, bitstring[i]);
        if(reversecodes.find(bit)!=reversecodes.end())
        {
           buffer.push_back(reversecodes[bit]);
           bit = "";
        }
    }

    sz = buffer.size();
    return buffer.data();
}


void writeHeader(char* path,map<unsigned char, string> codes,  int paddedBits){
    
    int size = codes.size();
    writeFileFromBuffer(path, (unsigned char*)&paddedBits, sizeof(int), 0);
    writeFileFromBuffer(path, (unsigned char*)&size, sizeof(int), 1);
    char nullBit = '\0';
    for(map<unsigned char, string>::iterator i = codes.begin(); i!=codes.end(); i++)
    {
        writeFileFromBuffer(path, (unsigned char*)&i->first, 1, 1);
        int len = i->second.size();
        writeFileFromBuffer(path, (unsigned char*)&len, sizeof(int), 1);
        writeFileFromBuffer(path, (unsigned char*)i->second.c_str(), i->second.size(), 1);
    }
}


unsigned char* readHeader(unsigned char* buffer, map<unsigned char, string> &codes, int& paddedBits, int &sz)
{
   paddedBits = *((int*)buffer);
   buffer = buffer+4;
   sz-=4;
   int size = *((int*)buffer);
   buffer = buffer+4;
   sz-=4;
   for(int i=0; i<size; i++)
   {    
       unsigned char key = buffer[0];
       buffer++;
       sz--;
       int len = *((int*)buffer);
       buffer+=4;
       sz-=4;
       char* value = (char*)malloc(len+1);

       for(int j = 0; j<len; j++)
       {
           value[j]=buffer[j];
       }
    //    value = (char*)buffer;
       buffer+=len;
       sz-=len;
       value[len]='\0';
       codes[key] = value;
   }

   return buffer;
}


void compressFile(char *path, std::vector<unsigned char> &data){
    int paddedBits = 0;
    vector<int> freqTable(256, 0);
    for (int i = 0; i < data.size(); i++)
    {
        freqTable[data[i]]++;
    }
    Tree *root = buildHuffmanTree(freqTable);
    map<unsigned char, string> codes;
    traverseHuffmanTree(root, "", "", codes);
    string outputString = getHuffmanBitstring(&data[0], codes, data.size(), paddedBits);
    
    int size  = outputString.size();
    vector<unsigned char> outputBufferV;
    getBufferFromString(outputString, outputBufferV, size);
    unsigned char* outputBuffer = outputBufferV.data();
    writeHeader(path, codes, paddedBits);
    writeFileFromBuffer(path, outputBuffer, size, 1);
}


void decompressFile( char* inputPath,  std::vector<unsigned char> &data)
{
    int sz = 0;
    map<unsigned char, string> codes;
    int paddedBits = 0;
    unsigned char* fileBuffer = readFileIntoBuffer(inputPath, sz);
    fileBuffer = readHeader(fileBuffer, codes, paddedBits, sz);
    string fileBitString = getStringFromBuffer(fileBuffer, sz);
    data.clear();
    getDecodedBuffer(fileBitString, data, codes, sz, paddedBits);
}