#include <iostream>
#include "gdal/gdal_priv.h"
#include "gdal/cpl_conv.h" // for CPLMalloc()
#include <string>
#include "main.h"
#include "EpiMap.h"

int main() {
    GDALAllRegister();
    GDALDriverManager *mgr = GetGDALDriverManager();
    GDALDriver *driverPng = mgr->GetDriverByName("PNG");
    GDALDataset *dataset;


    const char* filename = "N26E39.tif";
    dataset = (GDALDataset*) GDALOpen(filename, GA_Update);

    if (dataset != nullptr) {
        std::cout << "GDAL dataset opened!" << std::endl;
        int blockSizeX, blockSizeY;
        GDALRasterBand* band = dataset->GetRasterBand(1);
        band->GetBlockSize(&blockSizeX, &blockSizeY);
        printf("Block=%dx%d\n", blockSizeX, blockSizeY);
        uint8_t* buffer;
        int xSize = band->GetXSize();
        int ySize = band->GetYSize();
        printf("Band=%dx%d\n", xSize, ySize);
        GDALDataset* cropped = rasterPrepare(mgr, dataset, X_OFFSET, Y_OFFSET, WIDTH, HEIGHT);
        printf("Cropped to %dx%d\n", WIDTH, HEIGHT);

        printf("Loading data to the model...\n");
        auto model = new EpiMap(cropped->GetRasterBand(1));

        uint32_t n = model->getN();
        printf("Number of nodes: %d\n", n);
        printf("Building the k-d tree model...\n");

        model->buildSearchTree();
        point origin = {3000, 1000};
        double r = 50;
        std::vector<point> result = model->radiusSearch(origin, r);
        printf("Found %zu nodes!\n", result.size());

//        buffer = (uint8_t*) CPLMalloc(sizeof(uint8_t) * xSize * ySize);
//        rasterCopy(GF_Read, band, buffer, 6000, 6000, 2000, 4000);
//        std::cout << "Raster read!" << std::endl;
//        uint32_t loc = 4000*1200 + 123;
//        printf("Random pixel #%d: %d\n", loc, buffer[loc]);

//        std::cout << "Changing pixel..." << std::endl;
//        buffer[loc] = 23;
//        rasterCopy(GF_Write, band, buffer, 6000, 6000, 4000, 2000);
//        std::cout << "OK!" << std::endl;

        printf("Launching simulation...\n");
        model->simulate(1000, 0.05, 0.01, 0.01);

        GDALDataset* pngDs = driverPng->CreateCopy("./test.png", cropped, FALSE, nullptr, nullptr, nullptr);
        GDALFlushCache(pngDs);
        GDALClose(pngDs);
        GDALClose(cropped);
        GDALClose(dataset);

//        CPLFree(buffer);
        return 0;
    } else {
        std::cout << "Can't open dataset!" << std::endl;
        return 127;
    }
}

GDALColorTable *getColorTable() {
    auto colorTable = new GDALColorTable(GPI_RGB);
    GDALColorEntry col0, colS, colE, colI, colR;
    col0.c1 = 0, col0.c2 = 0, col0.c3 = 0, col0.c4 = 255;
    colS.c1 = 66, colS.c2 = 66, colS.c3 = 66, colS.c4 = 255;
    colE.c1 = 248, colE.c2 = 248, colE.c3 = 20, colE.c4 = 255;
    colI.c1 = 20, colI.c2 = 248, colI.c3 = 20, colI.c4 = 255;
    colR.c1 = 248, colR.c2 = 20, colR.c3 = 20, colR.c4 = 255;
    for (int i = 0; i < 1; i++) {
        colorTable->SetColorEntry(i, &col0);
    }
    for (int i = 1; i < 64; i++) {
        colorTable->SetColorEntry(i, &colS);
    }
    for (int i = 64; i < 128; i++) {
        colorTable->SetColorEntry(i, &colE);
    }
    for (int i = 128; i < 192; i++) {
        colorTable->SetColorEntry(i, &colI);
    }
    for (int i = 192; i < 256; i++) {
        colorTable->SetColorEntry(i, &colR);
    }
    return colorTable;
}

CPLErr applyColors(GDALDataset *dataset) {
    GDALColorTable* colorTable = getColorTable();
    return dataset->GetRasterBand(1)->SetColorTable(colorTable);
}

GDALDataset *
rasterPrepare(GDALDriverManager *mgr, GDALDataset *dataset, int xOffset, int yOffset, int width, int height) {
    GDALDriver *driverMem = mgr->GetDriverByName("MEM");
    GDALDataset* croppedRaster = driverMem->Create("", width, height, 1, GDT_Byte, nullptr);
    auto buffer = (uint8_t*) CPLMalloc(sizeof(uint8_t) * width * height);
    GDALRasterBand* inBand = dataset->GetRasterBand(1);
    GDALRasterBand* outBand = croppedRaster->GetRasterBand(1);
    rasterCopy(GF_Read, inBand, buffer, xOffset, yOffset, width, height);

    // filtering the buffer
    for (uint32_t i = 0; i < width * height; i++) {
        if (buffer[i] > 200) {
            buffer[i] = 1;
        } else {
            buffer[i] = 0;
        }
    }

    rasterCopy(GF_Write, outBand, buffer, 0, 0, width, height);
    applyColors(croppedRaster);
    CPLFree(buffer);
    return croppedRaster;
}
