/* Copyright (c) 2011 Khaled Mamou (kmamou at gmail dot com)
 All rights reserved.
 
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. The names of the contributors may not be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <list>
#include <cctype>
#include <cmath>

//#define _CRTDBG_MAP_ALLOC

#ifdef _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#include <stdlib.h>
#endif // _CRTDBG_MAP_ALLOC

#include "VHACD.h"
#include "oclHelper.h"

using namespace VHACD;
using namespace std;

class MyCallback : public IVHACD::IUserCallback {
public:
    MyCallback(void) {}
    ~MyCallback(){};
    void Update(const double overallProgress, const double stageProgress, const double operationProgress,
        const char* const stage, const char* const operation)
    {
        cout << setfill(' ') << setw(3) << (int)(overallProgress + 0.5) << "% "
             << "[ " << stage << " " << setfill(' ') << setw(3) << (int)(stageProgress + 0.5) << "% ] "
             << operation << " " << setfill(' ') << setw(3) << (int)(operationProgress + 0.5) << "%" << endl;
    };
};
class MyLogger : public IVHACD::IUserLogger {
public:
    MyLogger(void) {}
    MyLogger(const string& fileName) { OpenFile(fileName); }
    ~MyLogger(){};
    void Log(const char* const msg)
    {
        if (m_file.is_open()) {
            m_file << msg;
            m_file.flush();
        }
    }
    void OpenFile(const string& fileName)
    {
        m_file.open(fileName.c_str());
    }

private:
    ofstream m_file;
};
struct Material {

    float m_diffuseColor[3];
    float m_ambientIntensity;
    float m_specularColor[3];
    float m_emissiveColor[3];
    float m_shininess;
    float m_transparency;
    Material(void)
    {
        m_diffuseColor[0] = 0.5f;
        m_diffuseColor[1] = 0.5f;
        m_diffuseColor[2] = 0.5f;
        m_specularColor[0] = 0.5f;
        m_specularColor[1] = 0.5f;
        m_specularColor[2] = 0.5f;
        m_ambientIntensity = 0.4f;
        m_emissiveColor[0] = 0.0f;
        m_emissiveColor[1] = 0.0f;
        m_emissiveColor[2] = 0.0f;
        m_shininess = 0.4f;
        m_transparency = 0.5f;
    };
};
struct Parameters {
    unsigned int m_oclPlatformID;
    unsigned int m_oclDeviceID;
    string m_fileNameIn;
    string m_fileNameOut;
    string m_fileNameLog;
    bool m_run;
    IVHACD::Parameters m_paramsVHACD;
    Parameters(void)
    {
        m_run = true;
        m_oclPlatformID = 0;
        m_oclDeviceID = 0;
        m_fileNameIn = "";
        m_fileNameOut = "output.wrl";
        m_fileNameLog = "log.txt";
    }
};
bool LoadOFF(const string& fileName, vector<float>& points, vector<int>& triangles, IVHACD::IUserLogger& logger);
bool LoadOBJ(const string& fileName, vector<float>& points, vector<int>& triangles, IVHACD::IUserLogger& logger);
bool LoadSTL(const string& fileName, vector<float>& points, vector<int>& triangles, IVHACD::IUserLogger& logger);
bool SaveOFF(const string& fileName, const float* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, IVHACD::IUserLogger& logger);
bool SaveVRML2(ofstream& fout, const double* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, const Material& material, IVHACD::IUserLogger& logger);
bool SaveOBJ(ofstream& fout, const double* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, const Material& material, IVHACD::IUserLogger& logger, int convexPart, int vertexOffset);
bool SaveSTL(ofstream& fout, const double* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, IVHACD::IUserLogger& logger, int convexPart);
void GetFileExtension(const string& fileName, string& fileExtension);
void ComputeRandomColor(Material& mat);
void Usage(const Parameters& params);
void ParseParameters(int argc, char* argv[], Parameters& params);

#ifdef CL_VERSION_1_1
bool InitOCL(const unsigned int oclPlatformID, const unsigned int oclDeviceID, OCLHelper& oclHelper, std::ostringstream& msg);
#endif // CL_VERSION_1_1

int main(int argc, char* argv[])
{
    // --input camel.off --output camel_acd.wrl --log log.txt --resolution 1000000 --depth 20 --concavity 0.0025 --planeDownsampling 4 --convexhullDownsampling 4 --alpha 0.05 --beta 0.05 --gamma 0.00125 --pca 0 --mode 0 --maxNumVerticesPerCH 256 --minVolumePerCH 0.0001 --convexhullApproximation 1 --oclDeviceID 2
    {
        // set parameters
        Parameters params;
        ParseParameters(argc, argv, params);
        MyCallback myCallback;
        MyLogger myLogger(params.m_fileNameLog);
        params.m_paramsVHACD.m_logger = &myLogger;
        params.m_paramsVHACD.m_callback = &myCallback;
        Usage(params);
        if (!params.m_run) {
            return 0;
        }

        std::ostringstream msg;
#ifdef CL_VERSION_1_1
        msg << "+ OpenCL (ON)" << std::endl;
        OCLHelper oclHelper;
        if (params.m_paramsVHACD.m_oclAcceleration) {
            bool res = InitOCL(params.m_oclPlatformID,
                params.m_oclDeviceID,
                oclHelper,
                msg);
            if (!res) {
                myLogger.Log(msg.str().c_str());
                return -1;
            }
        }
#else //CL_VERSION_1_1
        msg << "+ OpenCL (OFF)" << std::endl;
#endif //CL_VERSION_1_1

#ifdef _OPENMP
        msg << "+ OpenMP (ON)" << std::endl;
#else
        msg << "+ OpenMP (OFF)" << std::endl;
#endif
        msg << "+ Parameters" << std::endl;
        msg << "\t input                                       " << params.m_fileNameIn << endl;
        msg << "\t resolution                                  " << params.m_paramsVHACD.m_resolution << endl;
        msg << "\t max. concavity                              " << params.m_paramsVHACD.m_concavity << endl;
        msg << "\t plane down-sampling                         " << params.m_paramsVHACD.m_planeDownsampling << endl;
        msg << "\t convex-hull down-sampling                   " << params.m_paramsVHACD.m_convexhullDownsampling << endl;
        msg << "\t alpha                                       " << params.m_paramsVHACD.m_alpha << endl;
        msg << "\t beta                                        " << params.m_paramsVHACD.m_beta << endl;
        msg << "\t maxhulls                                    " << params.m_paramsVHACD.m_maxConvexHulls << endl;
        msg << "\t pca                                         " << params.m_paramsVHACD.m_pca << endl;
        msg << "\t mode                                        " << params.m_paramsVHACD.m_mode << endl;
        msg << "\t max. vertices per convex-hull               " << params.m_paramsVHACD.m_maxNumVerticesPerCH << endl;
        msg << "\t min. volume to add vertices to convex-hulls " << params.m_paramsVHACD.m_minVolumePerCH << endl;
        msg << "\t convex-hull approximation                   " << params.m_paramsVHACD.m_convexhullApproximation << endl;
        msg << "\t OpenCL acceleration                         " << params.m_paramsVHACD.m_oclAcceleration << endl;
        msg << "\t OpenCL platform ID                          " << params.m_oclPlatformID << endl;
        msg << "\t OpenCL device ID                            " << params.m_oclDeviceID << endl;
        msg << "\t output                                      " << params.m_fileNameOut << endl;
        msg << "\t log                                         " << params.m_fileNameLog << endl;
        msg << "+ Load mesh" << std::endl;
        myLogger.Log(msg.str().c_str());

        cout << msg.str().c_str();

        // load mesh
        vector<float> points;
        vector<int> triangles;
        string fileExtension;
        GetFileExtension(params.m_fileNameIn, fileExtension);
        if (fileExtension == ".OFF") {
            if (!LoadOFF(params.m_fileNameIn, points, triangles, myLogger)) {
                return -1;
            }
        }
        else if (fileExtension == ".OBJ") {
            if (!LoadOBJ(params.m_fileNameIn, points, triangles, myLogger)) {
                return -1;
            }
        }
        else if (fileExtension == ".STL") {
            if (!LoadSTL(params.m_fileNameIn, points, triangles, myLogger)) {
                return -1;
            }
        }
        else {
            myLogger.Log("Format not supported!\n");
            return -1;
        }

        // run V-HACD
        IVHACD* interfaceVHACD = CreateVHACD();

#ifdef CL_VERSION_1_1
        if (params.m_paramsVHACD.m_oclAcceleration) {
            bool res = interfaceVHACD->OCLInit(oclHelper.GetDevice(), &myLogger);
            if (!res) {
                params.m_paramsVHACD.m_oclAcceleration = false;
            }
        }
#endif //CL_VERSION_1_1
        bool res = interfaceVHACD->Compute(&points[0], (unsigned int)points.size() / 3,
            (const uint32_t *)&triangles[0], (unsigned int)triangles.size() / 3, params.m_paramsVHACD);




        if (res) {
            std::string ext;
            if (params.m_fileNameOut.length() > 4) {
                ext = params.m_fileNameOut.substr(params.m_fileNameOut.length()-4);
            }
            std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
            if (ext != ".obj" && ext != ".stl") {
                // save output
                unsigned int nConvexHulls = interfaceVHACD->GetNConvexHulls();
                msg.str("");
                msg << "+ Generate output: " << nConvexHulls << " convex-hulls " << endl;
                myLogger.Log(msg.str().c_str());
                ofstream foutCH(params.m_fileNameOut.c_str());
                IVHACD::ConvexHull ch;
                if (foutCH.is_open()) {
                    Material mat;
                    for (unsigned int p = 0; p < nConvexHulls; ++p) {
                        interfaceVHACD->GetConvexHull(p, ch);
                        ComputeRandomColor(mat);
                        SaveVRML2(foutCH, ch.m_points, (const int *)ch.m_triangles, ch.m_nPoints, ch.m_nTriangles, mat, myLogger);
                        msg.str("");
                        msg << "\t CH[" << setfill('0') << setw(5) << p << "] " << ch.m_nPoints << " V, " << ch.m_nTriangles << " T" << endl;
                        myLogger.Log(msg.str().c_str());
                    }
                    foutCH.close();
                }
            }
            else if (ext == ".stl")
            {
                unsigned int nConvexHulls = interfaceVHACD->GetNConvexHulls();
                msg.str("");
                msg << "+ Generate output: " << nConvexHulls << " convex-hulls " << endl;
                myLogger.Log(msg.str().c_str());
                IVHACD::ConvexHull ch;
                if (params.m_fileNameOut.find('*') == string::npos)
                {
                    myLogger.Log("Could not find '*' in stl filename (did you try using quotation marks around the filename?)");
                    return -1;
                }
                for (unsigned int p = 0; p < nConvexHulls; ++p)
                {
                    string outname(params.m_fileNameOut);
                    char cp[100];
                    sprintf(cp, "%d", p);
                    outname.replace(outname.find('*'), 1, cp);
                    ofstream foutCH(outname.c_str());
                    if (foutCH.is_open())
                    {
                        interfaceVHACD->GetConvexHull(p, ch);
                        SaveSTL(foutCH, ch.m_points, (const int *)ch.m_triangles, ch.m_nPoints, ch.m_nTriangles, myLogger, p);
                        msg.str("");
                        msg << "\t CH[" << setfill('0') << setw(5) << p << "] " << ch.m_nPoints << " V, " << ch.m_nTriangles << " T" << endl;
                        myLogger.Log(msg.str().c_str());
                    foutCH.close();
                    }
                }
            }
            else {
                unsigned int nConvexHulls = interfaceVHACD->GetNConvexHulls();
                msg.str("");
                msg << "+ Generate output: " << nConvexHulls << " convex-hulls " << endl;
                myLogger.Log(msg.str().c_str());
                ofstream foutCH(params.m_fileNameOut.c_str());
                IVHACD::ConvexHull ch;
                if (foutCH.is_open()) {
                    Material mat;
                    int vertexOffset = 1;//obj wavefront starts counting at 1...
                    for (unsigned int p = 0; p < nConvexHulls; ++p) {
                        interfaceVHACD->GetConvexHull(p, ch);
                        SaveOBJ(foutCH, ch.m_points, (const int *)ch.m_triangles, ch.m_nPoints, ch.m_nTriangles, mat, myLogger, p, vertexOffset);
                        vertexOffset+=ch.m_nPoints;
                        msg.str("");
                        msg << "\t CH[" << setfill('0') << setw(5) << p << "] " << ch.m_nPoints << " V, " << ch.m_nTriangles << " T" << endl;
                        myLogger.Log(msg.str().c_str());
                    }
                    foutCH.close();
                }
            }
        }
        else {
            myLogger.Log("Decomposition cancelled by user!\n");
        }

#ifdef CL_VERSION_1_1
        if (params.m_paramsVHACD.m_oclAcceleration) {
            bool res = interfaceVHACD->OCLRelease(&myLogger);
            if (!res) {
                assert(-1);
            }
        }
#endif //CL_VERSION_1_1

        interfaceVHACD->Clean();
        interfaceVHACD->Release();
    }
#ifdef _CRTDBG_MAP_ALLOC
    _CrtDumpMemoryLeaks();
#endif // _CRTDBG_MAP_ALLOC
    return 0;
}

void Usage(const Parameters& params)
{
    std::ostringstream msg;
    msg << "V-HACD V" << VHACD_VERSION_MAJOR << "." << VHACD_VERSION_MINOR << endl;
    msg << "Syntax: testVHACD [options] --input infile.ext --output outfile.ext --log logfile.txt" << endl
        << endl;
    msg << "Options:" << endl;
    msg << "       --input                     Wavefront (.obj) or STL (.stl) input file name" << endl;
    msg << "       --output                    VRML 2.0 (.wrl) or multiple STL (*.stl, e.g. test_*.stl) output file name" << endl;
    msg << "       --log                       Log file name" << endl;
    msg << "       --resolution                Maximum number of voxels generated during the voxelization stage (default=100,000, range=10,000-16,000,000)" << endl;
    msg << "       --maxhulls                  Maximum number of convex hulls to produce." << endl;
    msg << "       --concavity                 Maximum allowed concavity (default=0.0025, range=0.0-1.0)" << endl;
    msg << "       --planeDownsampling         Controls the granularity of the search for the \"best\" clipping plane (default=4, range=1-16)" << endl;
    msg << "       --convexhullDownsampling    Controls the precision of the convex-hull generation process during the clipping plane selection stage (default=4, range=1-16)" << endl;
    msg << "       --alpha                     Controls the bias toward clipping along symmetry planes (default=0.05, range=0.0-1.0)" << endl;
    msg << "       --beta                      Controls the bias toward clipping along revolution axes (default=0.05, range=0.0-1.0)" << endl;
    msg << "       --gamma                     Controls the maximum allowed concavity during the merge stage (default=0.00125, range=0.0-1.0)" << endl;
    msg << "       --delta                     Controls the bias toward maximaxing local concavity (default=0.05, range=0.0-1.0)" << endl;
    msg << "       --pca                       Enable/disable normalizing the mesh before applying the convex decomposition (default=0, range={0,1})" << endl;
    msg << "       --mode                      0: voxel-based approximate convex decomposition, 1: tetrahedron-based approximate convex decomposition (default=0, range={0,1})" << endl;
    msg << "       --maxNumVerticesPerCH       Controls the maximum number of triangles per convex-hull (default=64, range=4-1024)" << endl;
    msg << "       --minVolumePerCH            Controls the adaptive sampling of the generated convex-hulls (default=0.0001, range=0.0-0.01)" << endl;
    msg << "       --convexhullApproximation   Enable/disable approximation when computing convex-hulls (default=1, range={0,1})" << endl;
    msg << "       --oclAcceleration           Enable/disable OpenCL acceleration (default=0, range={0,1})" << endl;
    msg << "       --oclPlatformID             OpenCL platform id (default=0, range=0-# OCL platforms)" << endl;
    msg << "       --oclDeviceID               OpenCL device id (default=0, range=0-# OCL devices)" << endl;
    msg << "       --help                      Print usage" << endl
        << endl;
    msg << "Examples:" << endl;
    msg << "       testVHACD.exe --input bunny.obj --output bunny_acd.wrl --log log.txt" << endl
        << endl;
    cout << msg.str();
    if (params.m_paramsVHACD.m_logger) {
        params.m_paramsVHACD.m_logger->Log(msg.str().c_str());
    }
}
void ParseParameters(int argc, char* argv[], Parameters& params)
{
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--input")) {
            if (++i < argc)
                params.m_fileNameIn = argv[i];
        }
        else if (!strcmp(argv[i], "--output")) {
            if (++i < argc)
                params.m_fileNameOut = argv[i];
        }
        else if (!strcmp(argv[i], "--log")) {
            if (++i < argc)
                params.m_fileNameLog = argv[i];
        }
        else if (!strcmp(argv[i], "--resolution")) {
            if (++i < argc)
                params.m_paramsVHACD.m_resolution = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--concavity")) {
            if (++i < argc)
                params.m_paramsVHACD.m_concavity = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--planeDownsampling")) {
            if (++i < argc)
                params.m_paramsVHACD.m_planeDownsampling = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--convexhullDownsampling")) {
            if (++i < argc)
                params.m_paramsVHACD.m_convexhullDownsampling = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--alpha")) {
            if (++i < argc)
                params.m_paramsVHACD.m_alpha = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--beta")) {
            if (++i < argc)
                params.m_paramsVHACD.m_beta = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--maxhulls")) {
            if (++i < argc)
                params.m_paramsVHACD.m_maxConvexHulls = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--pca")) {
            if (++i < argc)
                params.m_paramsVHACD.m_pca = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--mode")) {
            if (++i < argc)
                params.m_paramsVHACD.m_mode = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--maxNumVerticesPerCH")) {
            if (++i < argc)
                params.m_paramsVHACD.m_maxNumVerticesPerCH = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--minVolumePerCH")) {
            if (++i < argc)
                params.m_paramsVHACD.m_minVolumePerCH = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--convexhullApproximation")) {
            if (++i < argc)
                params.m_paramsVHACD.m_convexhullApproximation = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--oclAcceleration")) {
            if (++i < argc)
                params.m_paramsVHACD.m_oclAcceleration = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--oclPlatformID")) {
            if (++i < argc)
                params.m_oclPlatformID = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--oclDeviceID")) {
            if (++i < argc)
                params.m_oclDeviceID = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--help")) {
            params.m_run = false;
        }
    }
    params.m_paramsVHACD.m_resolution = (params.m_paramsVHACD.m_resolution < 64) ? 0 : params.m_paramsVHACD.m_resolution;
    params.m_paramsVHACD.m_planeDownsampling = (params.m_paramsVHACD.m_planeDownsampling < 1) ? 1 : params.m_paramsVHACD.m_planeDownsampling;
    params.m_paramsVHACD.m_convexhullDownsampling = (params.m_paramsVHACD.m_convexhullDownsampling < 1) ? 1 : params.m_paramsVHACD.m_convexhullDownsampling;
}

#ifdef CL_VERSION_1_1
bool InitOCL(const unsigned int oclPlatformID, const unsigned int oclDeviceID, OCLHelper& oclHelper, std::ostringstream& msg)
{

    bool res = true;
    vector<string> info;
    res = oclHelper.GetPlatformsInfo(info, "\t\t");
    if (!res)
        return res;

    const size_t numPlatforms = info.size();
    msg << "\t Number of OpenCL platforms: " << numPlatforms << endl;
    for (size_t i = 0; i < numPlatforms; ++i) {
        msg << "\t OpenCL platform [" << i << "]" << endl;
        msg << info[i];
    }
    msg << "\t Using OpenCL platform [" << oclPlatformID << "]" << endl;
    res = oclHelper.InitPlatform(oclPlatformID);
    if (!res)
        return res;

    info.clear();
    res = oclHelper.GetDevicesInfo(info, "\t\t");
    if (!res)
        return res;

    const size_t numDevices = info.size();
    msg << "\t Number of OpenCL devices: " << numDevices << endl;
    for (size_t i = 0; i < numDevices; ++i) {
        msg << "\t OpenCL device [" << i << "]" << endl;
        msg << info[i];
    }
    msg << "\t Using OpenCL device [" << oclDeviceID << "]" << endl;
    res = oclHelper.InitDevice(oclDeviceID);
    return res;
}
#endif // CL_VERSION_1_1
void GetFileExtension(const string& fileName, string& fileExtension)
{
    size_t lastDotPosition = fileName.find_last_of(".");
    if (lastDotPosition == string::npos) {
        fileExtension = "";
    }
    else {
        fileExtension = fileName.substr(lastDotPosition, fileName.size());
        transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::toupper);
    }
}
void ComputeRandomColor(Material& mat)
{
    mat.m_diffuseColor[0] = mat.m_diffuseColor[1] = mat.m_diffuseColor[2] = 0.0f;
    while (mat.m_diffuseColor[0] == mat.m_diffuseColor[1] || mat.m_diffuseColor[2] == mat.m_diffuseColor[1] || mat.m_diffuseColor[2] == mat.m_diffuseColor[0]) {
        mat.m_diffuseColor[0] = (rand() % 100) / 100.0f;
        mat.m_diffuseColor[1] = (rand() % 100) / 100.0f;
        mat.m_diffuseColor[2] = (rand() % 100) / 100.0f;
    }
}
bool LoadOFF(const string& fileName, vector<float>& points, vector<int>& triangles, IVHACD::IUserLogger& logger)
{
    FILE* fid = fopen(fileName.c_str(), "r");
    if (fid) {
        const string strOFF("OFF");
        char temp[1024];
        fscanf(fid, "%s", temp);
        if (string(temp) != strOFF) {
            logger.Log("Loading error: format not recognized \n");
            fclose(fid);
            return false;
        }
        else {
            int nv = 0;
            int nf = 0;
            int ne = 0;
            fscanf(fid, "%i", &nv);
            fscanf(fid, "%i", &nf);
            fscanf(fid, "%i", &ne);
            points.resize(nv * 3);
            triangles.resize(nf * 3);
            const int np = nv * 3;
            for (int p = 0; p < np; p++) {
                fscanf(fid, "%f", &(points[p]));
            }
            int s;
            for (int t = 0, r = 0; t < nf; ++t) {
                fscanf(fid, "%i", &s);
                if (s == 3) {
                    fscanf(fid, "%i", &(triangles[r++]));
                    fscanf(fid, "%i", &(triangles[r++]));
                    fscanf(fid, "%i", &(triangles[r++]));
                }
                else // Fix me: support only triangular meshes
                {
                    for (int h = 0; h < s; ++h)
                        fscanf(fid, "%i", &s);
                }
            }
            fclose(fid);
        }
    }
    else {
        logger.Log("Loading error: file not found \n");
        return false;
    }
    return true;
}
bool LoadOBJ(const string& fileName, vector<float>& points, vector<int>& triangles, IVHACD::IUserLogger& logger)
{
    const unsigned int BufferSize = 1024;
    FILE* fid = fopen(fileName.c_str(), "r");

    if (fid) {
        char buffer[BufferSize];
        int ip[4];
        float x[3];
        char* pch;
        char* str;
        while (!feof(fid)) {
            if (!fgets(buffer, BufferSize, fid)) {
                break;
            }
            else if (buffer[0] == 'v') {
                if (buffer[1] == ' ') {
                    str = buffer + 2;
                    for (int k = 0; k < 3; ++k) {
                        pch = strtok(str, " ");
                        if (pch)
                            x[k] = (float)atof(pch);
                        else {
                            return false;
                        }
                        str = NULL;
                    }
                    points.push_back(x[0]);
                    points.push_back(x[1]);
                    points.push_back(x[2]);
                }
            }
            else if (buffer[0] == 'f') {

                pch = str = buffer + 2;
                int k = 0;
                while (pch) {
                    pch = strtok(str, " ");
                    if (pch && *pch != '\n') {
                        ip[k++] = atoi(pch) - 1;
                    }
                    else {
                        break;
                    }
                    str = NULL;
                }
                if (k == 3) {
                    triangles.push_back(ip[0]);
                    triangles.push_back(ip[1]);
                    triangles.push_back(ip[2]);
                }
                else if (k == 4) {
                    triangles.push_back(ip[0]);
                    triangles.push_back(ip[1]);
                    triangles.push_back(ip[2]);

                    triangles.push_back(ip[0]);
                    triangles.push_back(ip[2]);
                    triangles.push_back(ip[3]);
                }
            }
        }
        fclose(fid);
    }
    else {
        logger.Log("File not found\n");
        return false;
    }
    return true;
}
char *nextword(char *str, char **next, IVHACD::IUserLogger& logger)
{
    char *start,*stop;

    start = &str[strspn(str," \t\n\v\f\r")];
    if (*start == '\0')
        return NULL;

    if (*start == '"' || *start == '\'')
    {
        stop = strchr(&start[1],*start);
        if (!stop)
        {
            logger.Log("Unbalanced quotes in input line\n");
            return NULL;
        }
        if (stop[1] && !isspace(stop[1]))
        {
            logger.Log("Input line quote not followed by whitespace\n");
            return NULL;
        }
        start++;
    }
    else
        stop = &start[strcspn(start," \t\n\v\f\r")];

    if (*stop == '\0')
        *next = NULL;
    else
        *next = stop+1;
    *stop = '\0';
    return start;
}
bool parse_line(char * line, char * copy, list<char *> &arg, unsigned int &narg, IVHACD::IUserLogger& logger)
{
    // duplicate line into copy string to break into words
    strcpy(copy,line);

    // strip any # comment by replacing it with 0
    // do not strip # inside single/double quotes

    char quote = '\0';
    char *ptr = copy;
    while (*ptr)
    {
        if (*ptr == '#' && !quote)
        {
            *ptr = '\0';
            break;
        }
        if (*ptr == quote)
            quote = '\0';
        else if (*ptr == '"' || *ptr == '\'')
            quote = *ptr;
        ptr++;
    }

    // point arg[] at each subsequent arg in copy string
    // nextword() inserts string terminators into copy string to delimit args
    // nextword() treats text between single/double quotes as one arg

    narg = 0;
    ptr = copy;
    char *next;

    while (ptr)
    {
        char * nw = nextword(ptr, &next, logger);
        if (!nw)
            return false;
        arg.push_back(nw);
        narg++;
        ptr = next;
    }
    return true;
}
bool meshtrifile_stl_binary(const string& fileName, vector<float>& points, vector<int>& triangles, IVHACD::IUserLogger& logger)
{
    unsigned int num_of_facets;
    std::ifstream stl_file;

    // open file for reading
    stl_file.open(fileName.c_str(), std::ifstream::in | std::ifstream::binary);

    // read 80 byte header into nirvana
    for (int i=0; i<20; i++){
      float dum;
      stl_file.read((char *)&dum, sizeof(float));
    }

    // read number of triangles
    stl_file.read((char *)&num_of_facets, sizeof(int));

    unsigned int nverts = 0;
    float tri_data[12];
    double minDist = -1.0;
    while(1) {
        // read one triangle dataset (normal + 3 vertices)
        for (int i=0; i<12; i++)
          stl_file.read((char *)&tri_data[i], sizeof(float));
        // read triangle attribute into nirvana
        short dum;
        stl_file.read((char *)&dum, sizeof(short));
        // error handling
        if (!stl_file)
        {
            logger.Log("Corrupt STL file: Error in reading binary STL file.\n");
            return false;
        }
        // add triangle to mesh
        // check if vertex is close to an existing one
        int vi[3] = {-1,-1,-1};
        for (int j = 0; j < 3; j++)
        {
            if (minDist < 0.)
                minDist = ((tri_data[3+0]-tri_data[6+0])*(tri_data[3+0]-tri_data[6+0])+(tri_data[3+1]-tri_data[6+1])*(tri_data[3+1]-tri_data[6+1])+(tri_data[3+2]-tri_data[6+2])*(tri_data[3+2]-tri_data[6+2]))*1e-10;

            bool found = false;
            for (unsigned int k =0; k < nverts*3; k+=3)
            {
                const double lenSqr = (tri_data[3*j+3+0]-points[k])*(tri_data[3*j+3+0]-points[k])+(tri_data[3*j+3+1]-points[k+1])*(tri_data[3*j+3+1]-points[k+1])+(tri_data[3*j+3+2]-points[k+2])*(tri_data[3*j+3+2]-points[k+2]);
                if (lenSqr < minDist)
                {
                    found = true;
                    vi[j] = k/3;
                    break;
                }
            }
            if(!found)
            {
                points.push_back(tri_data[3*j+3+0]);
                points.push_back(tri_data[3*j+3+1]);
                points.push_back(tri_data[3*j+3+2]);
                vi[j] = nverts;
                nverts++;
            }
        }
        triangles.push_back(vi[0]);
        triangles.push_back(vi[1]);
        triangles.push_back(vi[2]);

        if (triangles.size() == num_of_facets*3)
            break;
    }

    stl_file.close();
    return true;
}
bool LoadSTL(const string& fileName, vector<float>& points, vector<int>& triangles, IVHACD::IUserLogger& logger)
{
    const unsigned int BufferSize = 1024;
    FILE* fid = fopen(fileName.c_str(), "r");
    if (!fid)
    {
        logger.Log("File not found\n");
        return false;
    }

    unsigned int n, m;
    int iVertex = 0;
    double vertices[3][3];
    bool insideSolidObject = false;
    bool insideFacet = false;
    bool insideOuterLoop = false;
    double minDist = -1.0;

    int nLines = 0;
    char line[BufferSize], work[BufferSize];

    int nverts = 0;
    while (1)
    {
        // read a line from input script
        // n = length of line including str terminator, 0 if end of file
        // if line ends in continuation char '&', concatenate next line

        m = 0;
        while (1) {
            if (BufferSize-m < 2)
            {
                logger.Log("BufferSize too small\n");
                return false;
            }
            if (fgets(&line[m], BufferSize-m, fid) == NULL)
            {
                if (m)
                    n = strlen(line) + 1;
                else
                    n = 0;
                break;
            }
            m = strlen(line);
            if (line[m-1] != '\n')
                continue;

            m--;
            while (m >= 0 && isspace(line[m]))
                m--;
            if (m < 0 || line[m] != '&')
            {
                line[m+1] = '\0';
                n = m+2;
                break;
            }
        }

        // if n = 0, end-of-file
        // error if label_active is set, since label wasn't encountered
        // if original input file, code is done
        // else go back to previous input file

        if (n == 0)
            break;

        if (n > BufferSize)
        {
            logger.Log("BufferSize too small\n");
            return false;
        }

        // lines start with 1 (not 0)
        nLines++;

        unsigned int narg = 0;
        list<char *> arg;
        bool success = parse_line(line, work, arg, narg, logger);
        if (!success)
        {
            logger.Log("Could not parse STL file\n");
            return false;
        }

        // skip empty lines
        if(narg==0)
        {
            logger.Log("Note: Skipping empty line in STL file\n");
            continue;
        }

        char *firstWord = *(arg.begin());
        if (strcmp(firstWord,"solid") != 0 && nLines == 1)
        {
            fclose(fid);
            logger.Log("Note: solid keyword not found, assuming binary stl file\n");
            fid = NULL;
            return meshtrifile_stl_binary(fileName, points, triangles, logger);
        }

        // detect begin and end of a solid object, facet and vertices
        if (strcmp(firstWord,"solid") == 0)
        {
            if (insideSolidObject)
            {
                logger.Log("Corrupt or unknown STL file: New solid object begins without closing prior solid object.\n");
                return false;
            }
            insideSolidObject=true;
        }
        else if (strcmp(firstWord,"endsolid") == 0)
        {
            if (!insideSolidObject)
            {
                logger.Log("Corrupt or unknown STL file: End of solid object found, but no begin.\n");
                return false;
            }
            insideSolidObject=false;
        }
        // detect begin and end of a facet within a solids object
        else if (strcmp(firstWord,"facet") == 0)
        {
            if (insideFacet)
            {
                logger.Log("Corrupt or unknown STL file: New facet begins without closing prior facet.\n");
                return false;
            }
            if (!insideSolidObject)
            {
                logger.Log("Corrupt or unknown STL file: New facet begins outside solid object.\n");
                return false;
            }
            insideFacet = true;

            // check for keyword normal belonging to facet
            list<char*>::iterator it = arg.begin();
            it++;
            if (strcmp(*it,"normal") != 0)
            {
                logger.Log("Corrupt or unknown STL file: Facet normal not defined.\n");
                return false;
            }
            // do not import facet normal (is calculated later)
        }
        else if (strcmp(firstWord,"endfacet") == 0)
        {
            if (!insideFacet)
            {
                logger.Log("Corrupt or unknown STL file: End of facet found, but no begin.\n");
                return false;
            }
            insideFacet = false;
            if (iVertex != 3)
            {
                logger.Log("Corrupt or unknown STL file: Number of vertices not equal to three (no triangle).\n");
                return false;
            }

            // add triangle to mesh
            // check if vertex is close to an existing one
            int vi[3] = {-1,-1,-1};
            for (int j = 0; j < 3; j++)
            {
                if (minDist < 0.)
                {
                    minDist = ((vertices[0][0]-vertices[1][0])*(vertices[0][0]-vertices[1][0])+(vertices[0][1]-vertices[1][1])*(vertices[0][1]-vertices[1][1])+(vertices[0][2]-vertices[1][2])*(vertices[0][2]-vertices[1][2]))*1e-10;
                }

                bool found = false;
                for (int k =0; k < nverts*3; k+=3)
                {
                    const double lenSqr = (vertices[j][0]-points[k])*(vertices[j][0]-points[k])+(vertices[j][1]-points[k+1])*(vertices[j][1]-points[k+1])+(vertices[j][2]-points[k+2])*(vertices[j][2]-points[k+2]);
                    if (lenSqr < minDist)
                    {
                        found = true;
                        vi[j] = k/3;
                        break;
                    }
                }
                if(!found)
                {
                    points.push_back(vertices[j][0]);
                    points.push_back(vertices[j][1]);
                    points.push_back(vertices[j][2]);
                    vi[j] = nverts;
                    nverts++;
                }
            }
            triangles.push_back(vi[0]);
            triangles.push_back(vi[1]);
            triangles.push_back(vi[2]);
        }
        //detect begin and end of an outer loop within a facet
        else if (strcmp(firstWord,"outer") == 0)
        {
            if (insideOuterLoop)
            {
                logger.Log("Corrupt or unknown STL file: New outer loop begins without closing prior outer loop.\n");
                return false;
            }
            if (!insideFacet)
            {
                logger.Log("Corrupt or unknown STL file: New outer loop begins outside facet.\n");
                return false;
            }
            insideOuterLoop = true;
            iVertex = 0;
        }
        else if (strcmp(firstWord,"endloop") == 0)
        {
            if (!insideOuterLoop)
            {
                logger.Log("Corrupt or unknown STL file: End of outer loop found, but no begin.\n");
                return false;
            }
            insideOuterLoop=false;
        }

        else if (strcmp(firstWord,"vertex") == 0)
        {
            if (!insideOuterLoop)
            {
                logger.Log("Corrupt or unknown STL file: Vertex found outside a loop.\n");
                return false;
            }

            // read the vertex
            list<char*>::iterator it = arg.begin();
            for (int j=0; j<3; j++)
            {
                it++;
                vertices[iVertex][j]=atof(*it);
            }

            iVertex++;
            if (iVertex > 3)
            {
                logger.Log("Corrupt or unknown STL file: Can not have more than 3 vertices "
                           "in a facet (only triangular meshes supported).\n");
                return false;
            }
        }
    }
    return true;
}
bool SaveOFF(const string& fileName, const float* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, IVHACD::IUserLogger& logger)
{
    ofstream fout(fileName.c_str());
    if (fout.is_open()) {
        size_t nV = nPoints * 3;
        size_t nT = nTriangles * 3;
        fout << "OFF" << std::endl;
        fout << nPoints << " " << nTriangles << " " << 0 << std::endl;
        for (size_t v = 0; v < nV; v += 3) {
            fout << points[v + 0] << " "
                 << points[v + 1] << " "
                 << points[v + 2] << std::endl;
        }
        for (size_t f = 0; f < nT; f += 3) {
            fout << "3 " << triangles[f + 0] << " "
                 << triangles[f + 1] << " "
                 << triangles[f + 2] << std::endl;
        }
        fout.close();
        return true;
    }
    else {
        logger.Log("Can't open file\n");
        return false;
    }
}
bool SaveVRML2(ofstream& fout, const double* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, const Material& material, IVHACD::IUserLogger& logger)
{
    if (fout.is_open()) {
        fout.setf(std::ios::fixed, std::ios::floatfield);
        fout.setf(std::ios::showpoint);
        fout.precision(6);
        size_t nV = nPoints * 3;
        size_t nT = nTriangles * 3;
        fout << "#VRML V2.0 utf8" << std::endl;
        fout << "" << std::endl;
        fout << "# Vertices: " << nPoints << std::endl;
        fout << "# Triangles: " << nTriangles << std::endl;
        fout << "" << std::endl;
        fout << "Group {" << std::endl;
        fout << "    children [" << std::endl;
        fout << "        Shape {" << std::endl;
        fout << "            appearance Appearance {" << std::endl;
        fout << "                material Material {" << std::endl;
        fout << "                    diffuseColor " << material.m_diffuseColor[0] << " "
             << material.m_diffuseColor[1] << " "
             << material.m_diffuseColor[2] << std::endl;
        fout << "                    ambientIntensity " << material.m_ambientIntensity << std::endl;
        fout << "                    specularColor " << material.m_specularColor[0] << " "
             << material.m_specularColor[1] << " "
             << material.m_specularColor[2] << std::endl;
        fout << "                    emissiveColor " << material.m_emissiveColor[0] << " "
             << material.m_emissiveColor[1] << " "
             << material.m_emissiveColor[2] << std::endl;
        fout << "                    shininess " << material.m_shininess << std::endl;
        fout << "                    transparency " << material.m_transparency << std::endl;
        fout << "                }" << std::endl;
        fout << "            }" << std::endl;
        fout << "            geometry IndexedFaceSet {" << std::endl;
        fout << "                ccw TRUE" << std::endl;
        fout << "                solid TRUE" << std::endl;
        fout << "                convex TRUE" << std::endl;
        if (nV > 0) {
            fout << "                coord DEF co Coordinate {" << std::endl;
            fout << "                    point [" << std::endl;
            for (size_t v = 0; v < nV; v += 3) {
                fout << "                        " << points[v + 0] << " "
                     << points[v + 1] << " "
                     << points[v + 2] << "," << std::endl;
            }
            fout << "                    ]" << std::endl;
            fout << "                }" << std::endl;
        }
        if (nT > 0) {
            fout << "                coordIndex [ " << std::endl;
            for (size_t f = 0; f < nT; f += 3) {
                fout << "                        " << triangles[f + 0] << ", "
                     << triangles[f + 1] << ", "
                     << triangles[f + 2] << ", -1," << std::endl;
            }
            fout << "                ]" << std::endl;
        }
        fout << "            }" << std::endl;
        fout << "        }" << std::endl;
        fout << "    ]" << std::endl;
        fout << "}" << std::endl;
        return true;
    }
    else {
        logger.Log("Can't open file\n");
        return false;
    }
}


bool SaveOBJ(ofstream& fout, const double* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, const Material& material, IVHACD::IUserLogger& logger, int convexPart, int vertexOffset)
{
    if (fout.is_open()) {

        fout.setf(std::ios::fixed, std::ios::floatfield);
        fout.setf(std::ios::showpoint);
        fout.precision(6);
        size_t nV = nPoints * 3;
        size_t nT = nTriangles * 3;

		fout << "o convex_" << convexPart << std::endl;

        if (nV > 0) {
            for (size_t v = 0; v < nV; v += 3) {
                fout << "v " << points[v + 0] << " " << points[v + 1] << " " << points[v + 2] << std::endl;
            }
        }
        if (nT > 0) {
            for (size_t f = 0; f < nT; f += 3) {
                     fout << "f " 
                     << triangles[f + 0]+vertexOffset << " "
                     << triangles[f + 1]+vertexOffset << " "
                     << triangles[f + 2]+vertexOffset << " " << std::endl;
            }
        }
        return true;
    }
    else {
        logger.Log("Can't open file\n");
        return false;
    }
}

bool SaveSTL(ofstream& fout, const double* const& points, const int* const& triangles, const unsigned int& nPoints,
    const unsigned int& nTriangles, IVHACD::IUserLogger& logger, int convexPart)
{
    if (fout.is_open())
    {
        // write spaces int header
        char space = ' ';
        for (unsigned int i = 0; i < 80; i++)
            fout.write(&space, sizeof(char));

        fout.write((char *)&nTriangles, sizeof(int));

        for (unsigned int i = 0; i < nTriangles; i++)
        {
            float v[9];
            for (unsigned int j = 0; j < 3; j++)
            {
                v[3*j+0] = points[3*triangles[3*i+j]+0];
                v[3*j+1] = points[3*triangles[3*i+j]+1];
                v[3*j+2] = points[3*triangles[3*i+j]+2];
            }
            float n[3] = {0., 0., 0.};
            n[0] = (v[0+1]-v[6+1])*(v[3+2]-v[6+2])-(v[0+2]-v[6+2])*(v[3+1]-v[6+1]);
            n[1] = (v[0+2]-v[6+2])*(v[3+0]-v[6+0])-(v[0+0]-v[6+0])*(v[3+2]-v[6+2]);
            n[2] = (v[0+0]-v[6+0])*(v[3+1]-v[6+1])-(v[0+1]-v[6+1])*(v[3+0]-v[6+0]);
            const float invLen = 1.0/sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
            n[0] *= invLen;
            n[1] *= invLen;
            n[2] *= invLen;
            // write normal
            fout.write((char*)n, 3*sizeof(float));
            // write vertex data
            fout.write((char*)v, 9*sizeof(float));
            // write attribute
            short dum = 0;
            fout.write((char*)&dum, sizeof(short));
        }

        return true;
    }
    else {
        logger.Log("Can't open file\n");
        return false;
    }
}

