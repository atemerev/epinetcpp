//
// Created by atemerev on 5/8/20.
//

#ifndef EPINETCPP_MAIN_H
#define EPINETCPP_MAIN_H


#define WIDTH 4000
#define HEIGHT 2000
#define X_OFFSET 6000
#define Y_OFFSET 6000

CPLErr rasterCopy(GDALRWFlag direction, GDALRasterBand* band, uint8_t* buffer, int xOffset, int yOffset, int xSize, int ySize);

GDALColorTable* getColorTable();

CPLErr applyColors(GDALDataset* dataset);

GDALDataset* rasterPrepare(GDALDriverManager* mgr, GDALDataset* dataset, int xOffset, int yOffset, int width, int height);

#endif //EPINETCPP_MAIN_H
