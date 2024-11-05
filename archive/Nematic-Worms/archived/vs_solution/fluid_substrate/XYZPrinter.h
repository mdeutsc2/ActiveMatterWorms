
#ifndef __XYZ_PRINTER_H__
#define __XYZ_PRINTER_H__

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

//.. simulation interfacing
//#define _FRAMESPERFILE 120
const std::string key = "BS1_test";

//std::ifstream load;
//std::ofstream routput;
//const std::string xyzfile = "nwsim";
//const std::string xyz = ".xyz";
//const std::string xyzv = ".xyzv";
//const std::string xyzc = ".xyzc";
//const std::string csv = ".csv";
//const std::string dat = ".dat";
//const std::string datafolder = "./app_data/";
//-------------------------------------------------------------------------------------------------------------------

std::string CreateComment(float particlesPerWorm, float springConstantK2, float boxXSize, float boxYSize){
	stringstream line;
	line << particlesPerWorm << " " << springConstantK2 << " " << boxXSize << " " << boxYSize;
	return line.str();
}

// ------------------------------------------------------------------------------------------------------------------

class XYZPrinter{

	struct cxyz{
		char *c;
		float *x, *y, *z;
	};

	//.. output file stream
	std::ofstream fxyz;

	//.. comment line in every file
	std::string commentline;

	//. base file name
	std::string filename;

	//.. number of frames per file
	int framesperfile;

	//.. number of frame in current file
	int fileframecount;

	//.. current file number
	int filecount;

	//.. number of lines in queue
	int queuedlinecount;

	//.. queue state
	bool queued;

	//.. storage for data
	std::vector<char*> c;
	std::vector<float*> x, y, z;
	std::vector<int> sizes;

	//.. stores ids in allocated storage ptrs (c,x,y,z)
	std::vector<bool> storaged;
	
public:
	//.. not much here
	XYZPrinter();
	~XYZPrinter();

	//.. easy printing for user
	void Setup(const std::string &baseFileName, const std::string &commentLine, int printsPerFile);
	
	//.. update queue with user data which lives until print
	void AddToQueueReferences(int items, char *types, float *Rx, float *Ry, float *Rz);
	
	//.. update queue with user data to be temperarily stored in XYZPrinter instance
	void AddToQueueStorage(int items, char *types, float *Rx, float *Ry, float *Rz);

	//.. prints all queued data to file then dumbs references and storage
	void Print();

private:

	//.. interal helpers
	std::string ConstructFileName(const std::string &baseName, const std::string &ext);
	void OpenStream(const std::string &fileName);
	void CloseStream();
	void UpdateStream();
	void AddDataStorage(int items, char *t, float *x, float *y, float *z);
	void DumpAllData();
};

// -----------------------------------------------------------------------------------------------------------------------------------

XYZPrinter::XYZPrinter(){
	this->filecount = 0;
	this->fileframecount = 0;
	this->framesperfile = 1;
	this->queuedlinecount = 0;
	this->queued = false;
	this->commentline = "comment line";
	this->filename = "data";
}

XYZPrinter::~XYZPrinter(){
	if (queued) this->Print();
	if (fxyz.is_open()) this->CloseStream();
	this->DumpAllData();
}

void XYZPrinter::Setup(const std::string &fileName, const std::string &commentLine, int printsPerFile){
	this->filename = fileName;
	this->commentline = commentLine;
	this->framesperfile = printsPerFile;
	//this->fxyz.open(this->ConstructFileName(fileName, xyz), std::ios::out);
	this->OpenStream(this->ConstructFileName(fileName, xyz));
}

void XYZPrinter::AddToQueueReferences(int items = 0, char *types = NULL, float *Rx = NULL, float *Ry = NULL, float *Rz = NULL){
	if (items == 0) return;
	this->c.push_back(types);
	this->x.push_back(Rx);
	this->y.push_back(Ry);
	this->z.push_back(Rz);
	this->sizes.push_back(items);
	this->queuedlinecount += items;
	this->queued = true;
	this->storaged.push_back(false); // not in instance storage
}

void XYZPrinter::AddToQueueStorage(int items = 0, char *types = NULL, float *Rx = NULL, float *Ry = NULL, float *Rz = NULL){
	if (items == 0) return;
	this->AddDataStorage(items, types, Rx, Ry, Rz);
	this->sizes.push_back(items);
	this->queuedlinecount += items;
	this->queued = true;
	this->storaged.push_back(true);
}

void XYZPrinter::Print(){
	if (!this->queued) return;
	this->UpdateStream();

	fxyz << this->queuedlinecount << std::endl;
	fxyz << this->commentline << std::endl;
	for (int qid = 0; qid < sizes.size(); qid++){
		for (int lid = 0; lid < sizes[qid]; lid++){
			fxyz << this->c[qid][lid] << " "
				<< this->x[qid][lid] << " "
				<< this->y[qid][lid] << " "
				<< this->z[qid][lid] << std::endl;
		}
	}

	this->DumpAllData();
	this->queued = false;
}

// -----------------------------------------------------------------------------------------------------------------------------------

std::string XYZPrinter::ConstructFileName(const std::string &baseName, const std::string &ext){
	std::stringstream name;
	name << baseName << filecount++ << ext;
	return name.str();
}

void XYZPrinter::OpenStream(const std::string &fileName){
	if (!this->fxyz.is_open())
		this->fxyz.open(fileName, std::ios::out);
}

void XYZPrinter::CloseStream(){
	this->fxyz.close();
}

void XYZPrinter::UpdateStream(){
	if (this->fileframecount++ > this->framesperfile){
		this->CloseStream();
		this->OpenStream(this->ConstructFileName(this->filename, xyz));
	}
}

void XYZPrinter::AddDataStorage(int items, char *t, float *x, float *y, float *z){
	//.. allocate space
	this->c.push_back(new char[items]);
	this->x.push_back(new float[items]);
	this->y.push_back(new float[items]);
	this->z.push_back(new float[items]);

	//.. copy values with NULL safe assignment to 0.0f
	for (int i = 0; i < items; i++){
		if (t != NULL)	this->c.back()[i] = t[i];
		else this->c.back()[i] = 'x';
		if (x != NULL)	this->x.back()[i] = x[i];
		else this->x.back()[i] = 0.0f;
		if(y != NULL) this->y.back()[i] = y[i];
		else this->y.back()[i] = 0.0f;
		if(z != NULL) this->z.back()[i] = z[i];
		else this->z.back()[i] = 0.0f;
	}
}

void XYZPrinter::DumpAllData(){

	//.. free memory stored in instance
	for (int i = 0; i < this->sizes.size(); i++){
		if (this->storaged[i]){
			delete[] this->c[i];
			delete[] this->x[i];
			delete[] this->y[i];
			delete[] this->z[i];
		}
	}

	//.. clear data containers and flaggers
	this->c.clear();
	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->sizes.clear();
	this->storaged.clear();
}



#endif