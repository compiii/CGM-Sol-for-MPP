/**
 * @file compiii_4s.c
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-01-10
 *
 * @copyright Compiii Thesis Copyright (c) 2022
 *
 */
#include "compiii_4s.h"
#include "clogger.h"

int compiii4s()
{
    long d = 0;
    int i = 0, last = 0, a, b, iter;
    Block block;
    PixelatedBlock pixelatedBlock;
    Coord coord, coordBis;
    double buildTime;
    receivedList4s = createReceivedList4s();
    for (d = 1; (d <= maxDiag) && (i != maxEvalBlock); d++)
    {
        block = tabBlock[i];
        pixelatedBlock = tabPixelatedBlock[i];
        while (block.blockData.coreData.diag == d && (i != maxEvalBlock))
        {
            if (d == maxDiag)
                last = 1;

            if (d == 1)
            {
                for (a = pixelatedBlock.pBlockData.nbEvaluatePixel; a <= pixelatedBlock.pBlockData.maxDiag; a++)
                {
                    if (a <= pixelatedBlock.pBlockData.nbEvaluatePixel)
                    {
                        for (b = 1; b <= a; b++)
                        {
                            coord.i = a - b + 1;
                            coord.j = b;
                            if (d == 1000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                            {
                                logD("process %d start the pixel %d Coord(%d,%d) d = %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                            }
                            tmpDouble = MPI_Wtime();
                            computeMCOP(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], 0, 0, 0, 0);
                            calculTime += MPI_Wtime() - tmpDouble;
                            if (d == 1000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                            {
                                logD("process %d finish the pixel %d Coord(%d,%d) d = %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                            }

                            tmpDouble = MPI_Wtime();
                            BlockDataList *list1 = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                            list1 = deleteOccurBlockList(list1);
                            BlockDataList *list2 = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                            list2 = deleteOccurBlockList(list2);

                            BlockDataList *path = buildSendBlockDataList(list1, list2);

                            if (!isEmptyBlockDataList(path))
                                if (d == 10000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                                {
                                    logE("process %d send the pixel %d Coord(%d,%d) d = %d to process %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, path->blockData.rank, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                                }
                            if (a == pixelatedBlock.pBlockData.maxDiag)
                                sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 1);
                            else
                                sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 0);
                            buildTime += MPI_Wtime() - tmpDouble;
                        }
                    }
                    else
                    {
                        iter = 0;
                        for (b = 2 + (a - pixelatedBlock.pBlockData.nbEvaluatePixel - 1); b <= pixelatedBlock.pBlockData.nbEvaluatePixel; b++)
                        {
                            coord.i = pixelatedBlock.pBlockData.nbEvaluatePixel - iter++;
                            coord.j = b;
                            if (d == 1000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                            {
                                logD("process %d start the pixel %d Coord(%d,%d) d = %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                            }
                            tmpDouble = MPI_Wtime();
                            computeMCOP(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], 0, 0, 0, 1);
                            calculTime += MPI_Wtime() - tmpDouble;
                            if (d == 1000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                            {
                                logD("process %d finish the pixel %d Coord(%d,%d) d = %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                            }

                            tmpDouble = MPI_Wtime();
                            BlockDataList *list1 = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                            list1 = deleteOccurBlockList(list1);
                            BlockDataList *list2 = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                            list2 = deleteOccurBlockList(list2);

                            BlockDataList *path = buildSendBlockDataList(list1, list2);

                            if (!isEmptyBlockDataList(path))
                                if (d == 10000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                                {
                                    logE("process %d send the pixel %d Coord(%d,%d) d = %d to process %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, path->blockData.rank, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                                }
                            if (a == pixelatedBlock.pBlockData.maxDiag)
                                sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 1);
                            else
                                sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 0);
                            buildTime += MPI_Wtime() - tmpDouble;
                        }
                    }
                }
            }
            else if (d == 2)
            {
                for (a = 1; a <= pixelatedBlock.pBlockData.maxDiag; a++)
                {
                    if (a <= pixelatedBlock.pBlockData.nbEvaluatePixel)
                    {
                        for (b = 1; b <= a; b++)
                        {
                            coord.i = a - b + 1;
                            coord.j = b;

                            if (coord.i == 1)
                            {
                                PixelBlockDataList *path = pixelatedBlock.needPixelBlockDataListByColumn.pixelBlockDataList[coord.j];
                                if (!isEmptyPixelBlockDataList(path))
                                    while (path != NULL)
                                    {
                                        if (d == 2000 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                        {
                                            logW("process %d wait the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, path->pixelBlockData.coreData.id, path->pixelBlockData.coreData.coord.i, path->pixelBlockData.coreData.coord.j, tabBlockData[path->pixelBlockData.idBlockData].rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.coord.i, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.coord.j, d, MPI_Wtime() - startTime);
                                        }
                                        // tabs[pixelatedBlock.currentFrag].tabBlockData[path->pixelBlockData.idBlockData].fragLevel < maxFrag - 1
                                        if (tabBlockData[path->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[path->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].idBlockData].fragLevel)
                                        {
                                            coordBis.i = 1;
                                            coordBis.j = (int)ceil(coord.j / (float)pow(2, (tabBlockData[pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].idBlockData].fragLevel - tabBlockData[path->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                            PixelBlockData kBlock = tabPixelatedBlockData[path->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                            {
                                                receivePixelBlockData4s(kBlock, block.blockData);
                                            }
                                        }
                                        else
                                        {
                                            if (findTagToReceivedList4s(receivedList4s, path->pixelBlockData.coreData.id, path->pixelBlockData.coreData.address, path->pixelBlockData.coreData.coord) == 0)
                                            {
                                                receivePixelBlockData4s(path->pixelBlockData, block.blockData);
                                            }
                                        }
                                        if (d == 2000 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                        {
                                            logW("process %d receive the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, path->pixelBlockData.coreData.id, path->pixelBlockData.coreData.coord.i, path->pixelBlockData.coreData.coord.j, tabBlockData[path->pixelBlockData.idBlockData].rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.coord.i, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.coord.j, d, MPI_Wtime() - startTime);
                                        }
                                        path = path->next;
                                    }
                            }

                            if (coord.j == 1)
                            {
                                PixelBlockDataList *path = pixelatedBlock.needPixelBlockDataListByLine.pixelBlockDataList[coord.i];
                                if (!isEmptyPixelBlockDataList(path))
                                    while (path != NULL)
                                    {
                                        if (d == 2000 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                        {
                                            logW("process %d wait the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, path->pixelBlockData.coreData.id, path->pixelBlockData.coreData.coord.i, path->pixelBlockData.coreData.coord.j, tabBlockData[path->pixelBlockData.idBlockData].rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.coord.i, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.coord.j, d, MPI_Wtime() - startTime);
                                        }
                                        if (tabBlockData[path->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[path->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].idBlockData].fragLevel)
                                        {
                                            coordBis.i = (int)ceil(coord.i / (float)pow(2, (tabBlockData[pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].idBlockData].fragLevel - tabBlockData[path->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                            coordBis.j = 1;
                                            PixelBlockData kBlock = tabPixelatedBlockData[path->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                            {
                                                receivePixelBlockData4s(kBlock, block.blockData);
                                            }
                                        }
                                        else
                                        {
                                            if (findTagToReceivedList4s(receivedList4s, path->pixelBlockData.coreData.id, path->pixelBlockData.coreData.address, path->pixelBlockData.coreData.coord) == 0)
                                            {
                                                receivePixelBlockData4s(path->pixelBlockData, block.blockData);
                                            }
                                        }
                                        if (d == 2000 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                        {
                                            logW("process %d receive the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, path->pixelBlockData.coreData.id, path->pixelBlockData.coreData.coord.i, path->pixelBlockData.coreData.coord.j, tabBlockData[path->pixelBlockData.idBlockData].rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.coord.i, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.coord.j, d, MPI_Wtime() - startTime);
                                        }
                                        path = path->next;
                                    }
                            }

                            if (d == 2000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                            {
                                logD("process %d start the pixel %d Coord(%d,%d) d = %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                            }
                            tmpDouble = MPI_Wtime();
                            computeMCOP(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], 0, 0, 0, 0);
                            calculTime += MPI_Wtime() - tmpDouble;
                            if (d == 2000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                            {
                                logD("process %d finish the pixel %d Coord(%d,%d) d = %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                            }

                            tmpDouble = MPI_Wtime();
                            BlockDataList *list1 = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                            list1 = deleteOccurBlockList(list1);
                            BlockDataList *list2 = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                            list2 = deleteOccurBlockList(list2);

                            BlockDataList *path = buildSendBlockDataList(list1, list2);

                            if (!isEmptyBlockDataList(path))
                                if (d == 10000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                                {
                                    logE("process %d send the pixel %d Coord(%d,%d) d = %d to process %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, path->blockData.rank, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                                }
                            if (a == pixelatedBlock.pBlockData.maxDiag)
                                sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 1);
                            else
                                sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 0);
                            buildTime += MPI_Wtime() - tmpDouble;
                        }
                    }
                    else
                    {
                        iter = 0;
                        for (b = 2 + (a - pixelatedBlock.pBlockData.nbEvaluatePixel - 1); b <= pixelatedBlock.pBlockData.nbEvaluatePixel; b++)
                        {
                            coord.i = pixelatedBlock.pBlockData.nbEvaluatePixel - iter++;
                            coord.j = b;

                            if (d == 2000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                            {
                                logD("process %d start the pixel %d Coord(%d,%d) d = %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                            }
                            tmpDouble = MPI_Wtime();
                            computeMCOP(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], 0, 0, 0, 1);
                            calculTime += MPI_Wtime() - tmpDouble;

                            if (d == 2000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                            {
                                logD("process %d finish the pixel %d Coord(%d,%d) d = %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                            }

                            tmpDouble = MPI_Wtime();
                            BlockDataList *list1 = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                            list1 = deleteOccurBlockList(list1);
                            BlockDataList *list2 = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                            list2 = deleteOccurBlockList(list2);

                            BlockDataList *path = buildSendBlockDataList(list1, list2);

                            if (!isEmptyBlockDataList(path))
                                if (d == 10000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                                {
                                    logE("process %d send the pixel %d Coord(%d,%d) d = %d to process %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, path->blockData.rank, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                                }
                            if (a == pixelatedBlock.pBlockData.maxDiag)
                                sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 1);
                            else
                                sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 0);
                            buildTime += MPI_Wtime() - tmpDouble;
                        }
                    }
                }
            }
            else // if (1 == 2)
            {
                int finalize = 0, stop = 0, id1P, id2P, gap = 0, tmpGap = 0, diagEval = 1, step = 0, begin, end, init;
                PixelBlockData pixel, pixel1, pixel2;
                PixelBlockDataList *pathPbList;
                BlockDataList *pathBList;
                while (stop != 1)
                {
                    for (a = 1 + step; a <= pixelatedBlock.pBlockData.maxDiag && a <= diagEval; a++)
                    {
                        if (a <= pixelatedBlock.pBlockData.nbEvaluatePixel)
                        {
                            for (b = 1; b <= a; b++)
                            {
                                coord.i = a - b + 1;
                                coord.j = b;
                                pixel = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j];
                                id1P = id2P = ceil(pixel.coreData.id / 2.0);
                                if (pixel.coreData.id % 2 == 0)
                                {
                                    id1P += gap;
                                    id2P -= gap;
                                }
                                else if (pixel.coreData.id == 3)
                                {
                                    id2P -= 1;
                                }
                                else
                                {
                                    id1P += gap;
                                    id2P -= (gap + (gap != 0 ? 1 : 0));
                                }

                                if (id1P == pixel.coreData.id - 1)
                                    finalize = 1;

                                if (rank == 40000 && d == 4)
                                    logD("id1P = %d id2P = %d gap = %d value %d\n", id1P, id2P, gap, pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j]->pixelBlockData.coreData.id);

                                if (id1P <= pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j]->pixelBlockData.coreData.id)
                                {

                                    pathPbList = pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j];
                                    while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id1P)
                                    {
                                        pathPbList = pathPbList->next;
                                    }
                                    if (pathPbList == NULL)
                                    {
                                        logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id1P %d", rank, d, gap, id1P);
                                        exit(EXIT_FAILURE);
                                    }
                                    pixel1 = pathPbList->pixelBlockData;

                                    if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pixel1.idBlockData].rank)
                                        if (pixel1.idBlockData != block.blockData.coreData.id && (coord.i == 1)) //|| coord.j == 1
                                        {
                                            if (!isEmptyBlockDataList(block.needBlockDataList))
                                            {
                                                pathBList = block.needBlockDataList;
                                                while (pathBList != NULL && pixel1.idBlockData != pathBList->blockData.coreData.id)
                                                {
                                                    pathBList = pathBList->next;
                                                }
                                                if (pathBList != NULL)
                                                {
                                                    if (rank == 1000 && d == 6)
                                                    {
                                                        logD("id1P = %d id2P = %d gap = %d\n", id1P, id2P, gap);
                                                        printPixelBlockData(pixel1);
                                                    }

                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                    // checking before receiving
                                                    if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                    {
                                                        coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                        // printPixelBlockData(kBlock);
                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                        {
                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        receivePixelBlockData4s(pixel1, block.blockData);
                                                    }
                                                    // receivePixelBlockData(pixel1, block.blockData);
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d, of process %d %d) for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pixel1.coreData.id, pixel1.coreData.address, pixel1.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel1, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel1, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (findTagToReceivedList4s(receivedList4s, pixel1.coreData.id, pixel1.coreData.address, pixel1.coreData.coord) == 0)
                                                {
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                    // checking before receiving
                                                    if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                    {
                                                        coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                        {
                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        receivePixelBlockData4s(pixel1, block.blockData);
                                                    }
                                                    // receivePixelBlockData(pixel1, block.blockData);
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d, of process %d %d) for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                            }
                                        }

                                    if (pixel.coreData.id % 2 != 0 && id1P == id2P)
                                    {
                                        pathPbList = pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j];
                                        while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id1P - 1)
                                        {
                                            pathPbList = pathPbList->next;
                                        }
                                        if (pathPbList == NULL)
                                        {
                                            logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id1P - 1 = %d", rank, d, gap, id1P - 1);
                                            exit(EXIT_FAILURE);
                                        }

                                        if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pathPbList->pixelBlockData.idBlockData].rank)
                                            if (pathPbList->pixelBlockData.idBlockData != block.blockData.coreData.id && (coord.i == 1)) //|| coord.j == 1
                                            {
                                                if (!isEmptyBlockDataList(block.needBlockDataList))
                                                {
                                                    pathBList = block.needBlockDataList;
                                                    while (pathBList != NULL && pathPbList->pixelBlockData.idBlockData != pathBList->blockData.coreData.id)
                                                    {
                                                        pathBList = pathBList->next;
                                                    }
                                                    if (pathBList != NULL)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList4s(receivedList4s, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.address, pathPbList->pixelBlockData.coreData.coord) == 0)
                                                        {
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d wait the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                            // checking before receiving
                                                            if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                            {
                                                                coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                {
                                                                    receivePixelBlockData4s(kBlock, block.blockData);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                            }
                                                            // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d receive the pixel %d Coord(%d,%d) for  of process %d evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.address, pathPbList->pixelBlockData.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                }
                                else
                                {
                                    pixel1 = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i - (pixel.coreData.id - id1P)][coord.j];
                                }

                                /*if (id1P != id2P)
                                {*/
                                if (id2P <= pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i]->pixelBlockData.coreData.id)
                                {
                                    pathPbList = pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i];
                                    while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id2P)
                                    {
                                        if (d == 3000 && rank == 3)
                                        {
                                            printPixelBlockData(pathPbList->pixelBlockData);
                                        }
                                        pathPbList = pathPbList->next;
                                    }
                                    if (pathPbList == NULL)
                                    {
                                        logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id2P %d", rank, d, gap, id2P);
                                        exit(EXIT_FAILURE);
                                    }
                                    pixel2 = pathPbList->pixelBlockData;

                                    if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pixel2.idBlockData].rank)
                                        if (pixel2.idBlockData != block.blockData.coreData.id && (coord.j == 1)) // coord.i == 1 ||
                                        {
                                            if (!isEmptyBlockDataList(block.needBlockDataList))
                                            {
                                                pathBList = block.needBlockDataList;
                                                while (pathBList != NULL && pixel2.idBlockData != pathBList->blockData.coreData.id)
                                                {
                                                    if (rank == 200 && d == 4)
                                                    {
                                                        logD("id1P = %d id2P = %d pixel = %d pixel2 = %d gap = %d\n", id1P, id2P, pixel.coreData.id, pixel2.coreData.id, gap);
                                                        // printPixelBlockData(pathPbList->pixelBlockData);
                                                    }
                                                    if (rank == 4000 && d == 6)
                                                        printBlockData(pathBList->blockData);
                                                    pathBList = pathBList->next;
                                                }
                                                if (pathBList != NULL)
                                                {
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                    // checking before receiving
                                                    if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                    {

                                                        coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        if (rank == 70000 && d == 3 && pixel.coreData.id == 6)
                                                        {
                                                            printf("YEEEEEEEEEEEEEEEAH\n");
                                                            logD("id1P = %d id2P = %d gap = %d coord.i %d coord.j %d\n", id1P, id2P, gap, coordBis.i, coordBis.j);
                                                            printPixelBlockData(pixel2);
                                                        }
                                                        PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                        if (rank == 70000 && d == 3 && pixel.coreData.id == 6)
                                                        {

                                                            logD("id1P = %d id2P = %d gap = %d\n", id1P, id2P, gap);
                                                            printPixelBlockData(kBlock);
                                                            printf("YOOOOOOOOOOOOOOOOOAH\n");
                                                        }
                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                        {
                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        receivePixelBlockData4s(pixel2, block.blockData);
                                                    }
                                                    // receivePixelBlockData(pixel2, block.blockData);
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pixel2.coreData.id, pixel2.coreData.address, pixel2.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel2, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel2, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (findTagToReceivedList4s(receivedList4s, pixel2.coreData.id, pixel2.coreData.address, pixel2.coreData.coord) == 0)
                                                {
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                    // checking before receiving
                                                    if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                    {
                                                        coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                        {
                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        receivePixelBlockData4s(pixel2, block.blockData);
                                                    }
                                                    // receivePixelBlockData(pixel2, block.blockData);
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d, of process %d %d) for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                            }
                                        }

                                    if (pixel.coreData.id % 2 != 0 && id1P == id2P)
                                    {
                                        pathPbList = pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i];
                                        while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id2P - 1)
                                        {
                                            pathPbList = pathPbList->next;
                                        }
                                        if (pathPbList == NULL)
                                        {
                                            logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id2P - 1 = %d", rank, d, gap, id2P - 1);
                                            exit(EXIT_FAILURE);
                                        }

                                        if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pathPbList->pixelBlockData.idBlockData].rank)
                                            if (pathPbList->pixelBlockData.idBlockData != block.blockData.coreData.id && (coord.j == 1)) // coord.i == 1 ||
                                            {
                                                if (!isEmptyBlockDataList(block.needBlockDataList))
                                                {
                                                    pathBList = block.needBlockDataList;
                                                    while (pathBList != NULL && pathPbList->pixelBlockData.idBlockData != pathBList->blockData.coreData.id)
                                                    {
                                                        pathBList = pathBList->next;
                                                    }
                                                    if (pathBList != NULL)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList4s(receivedList4s, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.address, pathPbList->pixelBlockData.coreData.coord) == 0)
                                                        {
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d wait the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                            // checking before receiving
                                                            if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                            {
                                                                coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                {
                                                                    receivePixelBlockData4s(kBlock, block.blockData);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                            }
                                                            // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d receive the pixel %d Coord(%d,%d) for  of process %d evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.address, pathPbList->pixelBlockData.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                }
                                else
                                {
                                    pixel2 = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j - (pixel.coreData.id - id2P)];
                                }
                                /*}
                                else
                                {
                                    pixel2 = pixel1;
                                }*/

                                if (id1P == id2P)
                                {
                                    if (pixel.coreData.id % 2 == 0)
                                    {
                                        begin = pixel1.coreData.firstBound.i - 1;
                                        end = begin;
                                    }
                                    else
                                    {
                                        begin = pixel1.coreData.firstBound.i - 1;
                                        end = pixel1.coreData.firstBound.j;
                                    }
                                    init = 1;
                                    if (d == 4 && rank == 4000)
                                    {
                                        printPixelBlockData(pixel1);
                                        printPixelBlockData(pixel2);
                                        printf("begin %d end %d\n", begin, end);
                                    }
                                }
                                else
                                {
                                    if (pixel.coreData.id == 3)
                                    {
                                        begin = pixel1.coreData.firstBound.i - 1;
                                        end = begin;
                                        init = 1;
                                        /*printPixelBlockData(pixel1);
                                        printPixelBlockData(pixel2);
                                        printf("begin %d end %d\n", begin, end);*/
                                    }
                                    else
                                    {
                                        begin = pixel1.coreData.firstBound.i - 1;
                                        end = pixel1.coreData.firstBound.j - 1;
                                        init = 0;
                                        if (d == 50000)
                                        {
                                            printPixelBlockData(pixel1);
                                            printPixelBlockData(pixel2);
                                            printf("begin %d end %d id1P = %d id2P = %d gap = %d\n", begin, end, id1P, id2P, gap);
                                        }
                                    }
                                }

                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                {
                                    logD("process %d start the first (+,min) of the pixel %d Coord(%d,%d) (id1=%d; id2=%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, pixel1.coreData.id, pixel2.coreData.id, d, MPI_Wtime() - startTime);
                                }
                                tmpDouble = MPI_Wtime();
                                computeMCOP(pixel, begin, end, init, 2);
                                calculTime += MPI_Wtime() - tmpDouble;

                                if (d == 3 && rank == 7000000 && pixel.coreData.id == 6)
                                {
                                    printPixelBlockData(pixel1);
                                    printPixelBlockData(pixel2);
                                    printf("begin %d end %d\n", begin, end);
                                    FILE *file = fopen("test_p0_pixel_4s.csv", "a");
                                    for (int r = pixel.coreData.firstBound.i; r <= pixel.coreData.firstBound.j; r++)
                                    {
                                        for (int t = pixel.coreData.secondBound.i; t <= pixel.coreData.secondBound.j; t++)
                                        {
                                            // fprintf(file, "%d;", getMCOP(r, t));
                                            printf("%d\t", getMCOP(r, t));
                                        }
                                        // fprintf(file, "\n");
                                        printf("\n");
                                    }
                                    // fprintf(file, "\n");
                                    printf("\n");
                                    fclose(file);
                                }
                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                {
                                    logD("process %d finish the first (+,min) of the pixel %d Coord(%d,%d) (id1=%d; id2=%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, pixel1.coreData.id, pixel2.coreData.id, d, MPI_Wtime() - startTime);
                                }

                                if (id1P != id2P)
                                {
                                    if (id2P <= pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j]->pixelBlockData.coreData.id)
                                    {
                                        pathPbList = pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j];
                                        while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id2P)
                                        {
                                            pathPbList = pathPbList->next;
                                        }
                                        if (pathPbList == NULL)
                                        {
                                            logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id2P %d", rank, d, gap, id2P);
                                            exit(EXIT_FAILURE);
                                        }
                                        pixel1 = pathPbList->pixelBlockData;

                                        if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pixel1.idBlockData].rank)
                                            if (pixel1.idBlockData != block.blockData.coreData.id && (coord.i == 1)) //|| coord.j == 1
                                            {
                                                if (!isEmptyBlockDataList(block.needBlockDataList))
                                                {
                                                    pathBList = block.needBlockDataList;
                                                    while (pathBList != NULL && pixel1.idBlockData != pathBList->blockData.coreData.id)
                                                    {
                                                        pathBList = pathBList->next;
                                                    }
                                                    if (pathBList != NULL)
                                                    {

                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel1, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel1, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        if (d == 8000 && rank == 0 && pixel1.coreData.id == 2)
                                                        {
                                                            printf("HELOOOOOOOOOOO\nRECIEVE LIST BE\n");
                                                            ReceivedList4s *path = receivedList4s;
                                                            while (path != NULL)
                                                            {
                                                                printf("(id %d, tag %d) \t", path->id, path->tag);
                                                                path = path->next;
                                                            }
                                                            printf(" - RECIEVE LIST END\n");
                                                        }
                                                        if (findTagToReceivedList4s(receivedList4s, pixel1.coreData.id, pixel1.coreData.address, pixel1.coreData.coord) == 0)
                                                        {
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d wait the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                            // checking before receiving
                                                            if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                            {
                                                                coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                {
                                                                    receivePixelBlockData4s(kBlock, block.blockData);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                receivePixelBlockData4s(pixel1, block.blockData);
                                                            }
                                                            // receivePixelBlockData(pixel1, block.blockData);
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d receive the pixel %d Coord(%d,%d) for  of process %d evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pixel1.coreData.id, pixel1.coreData.address, pixel1.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel1, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel1, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                    else
                                    {
                                        pixel1 = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i - (pixel.coreData.id - id2P)][coord.j];
                                    }

                                    if (rank == 300 && d == 3)
                                    {
                                        logD("id1P = %d id2P = %d pixel = %d need pixel = %d gap = %d\n", id1P, id2P, pixel.coreData.id, pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i]->pixelBlockData.coreData.id, gap);
                                    }
                                    if (id1P <= pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i]->pixelBlockData.coreData.id)
                                    {
                                        pathPbList = pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i];
                                        while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id1P)
                                        {
                                            pathPbList = pathPbList->next;
                                        }
                                        if (pathPbList == NULL)
                                        {
                                            logE("the list of need pixel of blocks must contains the id %d", id1P);
                                            exit(EXIT_FAILURE);
                                        }
                                        pixel2 = pathPbList->pixelBlockData;

                                        if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pixel2.idBlockData].rank)
                                            if (pixel2.idBlockData != block.blockData.coreData.id && (coord.j == 1)) // coord.i == 1 ||
                                            {
                                                if (!isEmptyBlockDataList(block.needBlockDataList))
                                                {
                                                    pathBList = block.needBlockDataList;
                                                    while (pathBList != NULL && pixel2.idBlockData != pathBList->blockData.coreData.id)
                                                    {
                                                        pathBList = pathBList->next;
                                                    }
                                                    if (pathBList != NULL)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            if (rank == 4000 && d == 4)
                                                            {
                                                                logD("id1P = %d id2P = %d gap = %d\n", id1P, id2P, gap);
                                                                printPixelBlockData(pixel2);
                                                                printBlockData(tabBlockData[pixel2.idBlockData]);
                                                            }
                                                            receivePixelBlockData4s(pixel2, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel2, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList4s(receivedList4s, pixel2.coreData.id, pixel2.coreData.address, pixel2.coreData.coord) == 0)
                                                        {
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d wait the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                            // checking before receiving
                                                            if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                            {
                                                                coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                {
                                                                    receivePixelBlockData4s(kBlock, block.blockData);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                receivePixelBlockData4s(pixel2, block.blockData);
                                                            }
                                                            // receivePixelBlockData(pixel2, block.blockData);
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d receive the pixel %d Coord(%d,%d) for  of process %d evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pixel2.coreData.id, pixel2.coreData.address, pixel2.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel2, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel2, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                    else
                                    {
                                        pixel2 = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j - (pixel.coreData.id - id1P)];
                                    }

                                    begin = pixel2.coreData.secondBound.i;
                                    end = pixel2.coreData.secondBound.j;
                                    init = 0;
                                    if (d == 50000)
                                    {
                                        printPixelBlockData(pixel1);
                                        printPixelBlockData(pixel2);
                                        printf("begin %d end %d\n", begin, end);
                                    }
                                    /*printPixelBlockData(pixel1);
                                    printPixelBlockData(pixel2);
                                    printf("begin %d end %d\n", begin, end);*/

                                    if (d == 8 && verboseDebug && rank == 10000) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d start the second (+,min) of the pixel %d Coord(%d,%d) (id1=%d; id2=%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, pixel1.coreData.id, pixel2.coreData.id, d, MPI_Wtime() - startTime);
                                    }
                                    tmpDouble = MPI_Wtime();
                                    computeMCOP(pixel, begin, end, init, 2);
                                    calculTime += MPI_Wtime() - tmpDouble;
                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d finish the second (+,min) of the pixel %d Coord(%d,%d) (id1=%d; id2=%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, pixel1.coreData.id, pixel2.coreData.id, d, MPI_Wtime() - startTime);
                                    }
                                }

                                if (d == 3 && rank == 7000000 && pixel.coreData.id == 6)
                                {
                                    printPixelBlockData(pixel1);
                                    printPixelBlockData(pixel2);
                                    if (pixel1.coreData.id == 2 && pixel2.coreData.id == 12)
                                    {
                                        for (int r = pixel1.coreData.firstBound.i; r <= pixel1.coreData.firstBound.j; r++)
                                        {
                                            for (int t = pixel1.coreData.secondBound.i; t <= pixel1.coreData.secondBound.j; t++)
                                            {
                                                // fprintf(file, "%d;", getMCOP(r, t));
                                                printf("%d\t", getMCOP(r, t));
                                            }
                                            // fprintf(file, "\n");
                                            printf("\n");
                                        }
                                        // fprintf(file, "\n");
                                        printf("\n");
                                        for (int r = pixel2.coreData.firstBound.i; r <= pixel2.coreData.firstBound.j; r++)
                                        {
                                            for (int t = pixel2.coreData.secondBound.i; t <= pixel2.coreData.secondBound.j; t++)
                                            {
                                                // fprintf(file, "%d;", getMCOP(r, t));
                                                printf("%d\t", getMCOP(r, t));
                                            }
                                            // fprintf(file, "\n");
                                            printf("\n");
                                        }
                                        // fprintf(file, "\n");
                                        printf("\n");
                                    }
                                    printf("begin %d end %d\n", begin, end);
                                    FILE *file = fopen("test_p0_pixel_4s.csv", "a");
                                    for (int r = pixel.coreData.firstBound.i; r <= pixel.coreData.firstBound.j; r++)
                                    {
                                        for (int t = pixel.coreData.secondBound.i; t <= pixel.coreData.secondBound.j; t++)
                                        {
                                            // fprintf(file, "%d;", getMCOP(r, t));
                                            printf("%d\t", getMCOP(r, t));
                                        }
                                        // fprintf(file, "\n");
                                        printf("\n");
                                    }
                                    // fprintf(file, "\n");
                                    printf("\n");
                                    fclose(file);
                                }
                                if (finalize)
                                {
                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d start the finalisation of the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                    }

                                    /*if (d == 8 && rank == 4000)
                                    {
                                        printPixelBlockData(pixel);
                                        int r, t;
                                        for (r = pixel.coreData.firstBound.i; r <= pixel.coreData.firstBound.j; r++)
                                        {
                                            for (t = pixel.coreData.secondBound.i; t <= pixel.coreData.secondBound.j; t++)
                                            {
                                                printf("%d\t", getMCOP(r, t));
                                            }
                                            printf("\n");
                                        }
                                    }*/

                                    tmpDouble = MPI_Wtime();
                                    computeMCOP(pixel, 0, 0, 0, 3);
                                    calculTime += MPI_Wtime() - tmpDouble;

                                    if (d == 3 && rank == 7000000 && pixel.coreData.id == 6)
                                    {
                                        printPixelBlockData(pixel1);
                                        printPixelBlockData(pixel2);
                                        printf("begin %d end %d\n", begin, end);
                                        FILE *file = fopen("test_p0_pixel_4s.csv", "a");
                                        for (int r = pixel.coreData.firstBound.i; r <= pixel.coreData.firstBound.j; r++)
                                        {
                                            for (int t = pixel.coreData.secondBound.i; t <= pixel.coreData.secondBound.j; t++)
                                            {
                                                // fprintf(file, "%d;", getMCOP(r, t));
                                                printf("%d\t", getMCOP(r, t));
                                            }
                                            // fprintf(file, "\n");
                                            printf("\n");
                                        }
                                        // fprintf(file, "\n");
                                        printf("\n");
                                        fclose(file);
                                    }

                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d finishs the finalisation of the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                    }

                                    tmpDouble = MPI_Wtime();
                                    BlockDataList *list1 = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                                    list1 = deleteOccurBlockList(list1);
                                    BlockDataList *list2 = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                                    list2 = deleteOccurBlockList(list2);

                                    BlockDataList *path = buildSendBlockDataList(list1, list2);

                                    if (!isEmptyBlockDataList(path))
                                        if (d == 10000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                                        {
                                            logE("process %d send the pixel %d Coord(%d,%d) d = %d to process %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, path->blockData.rank, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                                        }
                                    if (a == pixelatedBlock.pBlockData.maxDiag)
                                        sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 1);
                                    else
                                        sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 0);
                                    buildTime += MPI_Wtime() - tmpDouble;
                                }
                            }
                        }
                        else
                        {
                            iter = 0;
                            for (b = 2 + (a - pixelatedBlock.pBlockData.nbEvaluatePixel - 1); b <= pixelatedBlock.pBlockData.nbEvaluatePixel; b++)
                            {
                                coord.i = pixelatedBlock.pBlockData.nbEvaluatePixel - iter++;
                                coord.j = b;

                                pixel = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j];
                                id1P = id2P = ceil(pixel.coreData.id / 2.0);
                                if (pixel.coreData.id % 2 == 0)
                                {
                                    id1P += gap;
                                    id2P -= gap;
                                }
                                else
                                {
                                    id1P += gap;
                                    id2P -= (gap + (gap != 0 ? 1 : 0));
                                }

                                if (id1P == pixel.coreData.id - 1)
                                    finalize = 1;

                                if (id1P <= pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j]->pixelBlockData.coreData.id)
                                {
                                    pathPbList = pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j];
                                    while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id1P)
                                    {
                                        pathPbList = pathPbList->next;
                                    }
                                    if (pathPbList == NULL)
                                    {
                                        logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id1P %d", rank, d, gap, id1P);
                                        exit(EXIT_FAILURE);
                                    }
                                    pixel1 = pathPbList->pixelBlockData;

                                    if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pixel1.idBlockData].rank)
                                        if (pixel1.idBlockData != block.blockData.coreData.id && (coord.i == 1)) //|| coord.j == 1
                                        {
                                            if (!isEmptyBlockDataList(block.needBlockDataList))
                                            {
                                                pathBList = block.needBlockDataList;
                                                while (pathBList != NULL && pixel1.idBlockData != pathBList->blockData.coreData.id)
                                                {
                                                    pathBList = pathBList->next;
                                                }
                                                if (pathBList != NULL)
                                                {
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }

                                                    // checking before receiving
                                                    if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                    {
                                                        coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                        {
                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        receivePixelBlockData4s(pixel1, block.blockData);
                                                    }
                                                    // receivePixelBlockData(pixel1, block.blockData);
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d, of process %d %d) for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pixel1.coreData.id, pixel1.coreData.address, pixel1.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }

                                                        // checking before receiving
                                                        if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel1, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel1, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (findTagToReceivedList4s(receivedList4s, pixel1.coreData.id, pixel1.coreData.address, pixel1.coreData.coord) == 0)
                                                {
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }

                                                    // checking before receiving
                                                    if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                    {
                                                        coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                        {
                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        receivePixelBlockData4s(pixel1, block.blockData);
                                                    }
                                                    // receivePixelBlockData(pixel1, block.blockData);
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d, of process %d %d) for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                            }

                                            if (pixel.coreData.id % 2 != 0 && id1P == id2P)
                                            {
                                                pathPbList = pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j];
                                                while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id1P - 1)
                                                {
                                                    pathPbList = pathPbList->next;
                                                }
                                                if (pathPbList == NULL)
                                                {
                                                    logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id1P - 1 = %d", rank, d, gap, id1P - 1);
                                                    exit(EXIT_FAILURE);
                                                }

                                                if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pathPbList->pixelBlockData.idBlockData].rank)
                                                    if (pathPbList->pixelBlockData.idBlockData != block.blockData.coreData.id && (coord.i == 1)) //|| coord.j == 1
                                                    {
                                                        if (!isEmptyBlockDataList(block.needBlockDataList))
                                                        {
                                                            pathBList = block.needBlockDataList;
                                                            while (pathBList != NULL && pathPbList->pixelBlockData.idBlockData != pathBList->blockData.coreData.id)
                                                            {
                                                                pathBList = pathBList->next;
                                                            }
                                                            if (pathBList != NULL)
                                                            {
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d wait the pixel %d Coord(%d,%d) for ev of process %d aluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }

                                                                // checking before receiving
                                                                if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                                {
                                                                    coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                    coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                    PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                    if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                    {
                                                                        receivePixelBlockData4s(kBlock, block.blockData);
                                                                    }
                                                                }
                                                                else
                                                                {
                                                                    receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                                }
                                                                // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d receive the pixel %d Coord(%d,%d) for eval of process %d uate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                if (findTagToReceivedList4s(receivedList4s, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.address, pathPbList->pixelBlockData.coreData.coord) == 0)
                                                                {
                                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                    {
                                                                        logW("process %d wait the pixel %d Coord(%d,%d) for evaluate of process %d  the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                    }

                                                                    // checking before receiving
                                                                    if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                                    {
                                                                        coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                        coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                        PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                        {
                                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                                        }
                                                                    }
                                                                    else
                                                                    {
                                                                        receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                                    }
                                                                    // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                    {
                                                                        logW("process %d receive the pixel %d Coord(%d,%d) for evaluate of process %d  the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                            if (findTagToReceivedList4s(receivedList4s, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.address, pathPbList->pixelBlockData.coreData.coord) == 0)
                                                            {
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d wait the pixel %d Coord(%d,%d) for ev of process %d aluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }

                                                                // checking before receiving
                                                                if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                                {
                                                                    coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                    coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                    PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                    if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                    {
                                                                        receivePixelBlockData4s(kBlock, block.blockData);
                                                                    }
                                                                }
                                                                else
                                                                {
                                                                    receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                                }
                                                                // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d receive the pixel %d Coord(%d,%d) for eval of process %d uate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }
                                                            }
                                                        }
                                                    }
                                            }
                                        }
                                }
                                else
                                {
                                    pixel1 = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i - (pixel.coreData.id - id1P)][coord.j];
                                }

                                if (id2P <= pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i]->pixelBlockData.coreData.id)
                                {
                                    pathPbList = pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i];
                                    while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id2P)
                                    {
                                        pathPbList = pathPbList->next;
                                    }
                                    if (pathPbList == NULL)
                                    {
                                        logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id2P %d", rank, d, gap, id2P);
                                        exit(EXIT_FAILURE);
                                    }
                                    pixel2 = pathPbList->pixelBlockData;

                                    if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pixel2.idBlockData].rank)
                                        if (pixel2.idBlockData != block.blockData.coreData.id && (coord.j == 1)) // coord.i == 1 ||
                                        {
                                            if (!isEmptyBlockDataList(block.needBlockDataList))
                                            {
                                                pathBList = block.needBlockDataList;
                                                while (pathBList != NULL && pixel2.idBlockData != pathBList->blockData.coreData.id)
                                                {
                                                    pathBList = pathBList->next;
                                                }
                                                if (pathBList != NULL)
                                                {
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }

                                                    // checking before receiving
                                                    if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                    {
                                                        coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                        {
                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        receivePixelBlockData4s(pixel2, block.blockData);
                                                    }
                                                    // receivePixelBlockData(pixel2, block.blockData);
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d, of process %d %d) for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pixel2.coreData.id, pixel2.coreData.address, pixel2.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel2, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel2, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (findTagToReceivedList4s(receivedList4s, pixel2.coreData.id, pixel2.coreData.address, pixel2.coreData.coord) == 0)
                                                {
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                    // checking before receiving
                                                    if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                    {
                                                        coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                        PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                        {
                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        receivePixelBlockData4s(pixel2, block.blockData);
                                                    }
                                                    // receivePixelBlockData(pixel2, block.blockData);
                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d, of process %d %d) for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                            }

                                            if (pixel.coreData.id % 2 != 0 && id1P == id2P)
                                            {
                                                pathPbList = pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i];
                                                while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id2P - 1)
                                                {
                                                    pathPbList = pathPbList->next;
                                                }
                                                if (pathPbList == NULL)
                                                {
                                                    logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id2P - 1 = %d", rank, d, gap, id2P - 1);
                                                    exit(EXIT_FAILURE);
                                                }

                                                if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pathPbList->pixelBlockData.idBlockData].rank)
                                                    if (pathPbList->pixelBlockData.idBlockData != block.blockData.coreData.id && (coord.j == 1)) // coord.i == 1 ||
                                                    {
                                                        if (!isEmptyBlockDataList(block.needBlockDataList))
                                                        {
                                                            pathBList = block.needBlockDataList;
                                                            while (pathBList != NULL && pathPbList->pixelBlockData.idBlockData != pathBList->blockData.coreData.id)
                                                            {
                                                                pathBList = pathBList->next;
                                                            }
                                                            if (pathBList != NULL)
                                                            {
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d wait the pixel %d Coord(%d,%d) for ev of process %d aluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }
                                                                // checking before receiving
                                                                if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                                {
                                                                    coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                    coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                    PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                    if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                    {
                                                                        receivePixelBlockData4s(kBlock, block.blockData);
                                                                    }
                                                                }
                                                                else
                                                                {
                                                                    receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                                }
                                                                // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d wait the pixel %d Coord(%d,%d) for ev of process %d aluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d receive the pixel %d Coord(%d,%d) for eval of process %d uate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                if (findTagToReceivedList4s(receivedList4s, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.address, pathPbList->pixelBlockData.coreData.coord) == 0)
                                                                {
                                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                    {
                                                                        logW("process %d wait the pixel %d Coord(%d,%d) for evaluate of process %d  the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                    }
                                                                    // checking before receiving
                                                                    if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                                    {
                                                                        coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                        coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                        PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                        if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                        {
                                                                            receivePixelBlockData4s(kBlock, block.blockData);
                                                                        }
                                                                    }
                                                                    else
                                                                    {
                                                                        receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                                    }
                                                                    // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                    {
                                                                        logW("process %d receive the pixel %d Coord(%d,%d) for evaluate of process %d  the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                            if (findTagToReceivedList4s(receivedList4s, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.address, pathPbList->pixelBlockData.coreData.coord) == 0)
                                                            {
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d wait the pixel %d Coord(%d,%d) for ev of process %d aluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }
                                                                // checking before receiving
                                                                if (tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                                {
                                                                    coordBis.i = (int)ceil(pathPbList->pixelBlockData.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                    coordBis.j = (int)ceil(pathPbList->pixelBlockData.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pathPbList->pixelBlockData.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                    PixelBlockData kBlock = tabPixelatedBlockData[pathPbList->pixelBlockData.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                    if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                    {
                                                                        receivePixelBlockData4s(kBlock, block.blockData);
                                                                    }
                                                                }
                                                                else
                                                                {
                                                                    receivePixelBlockData4s(pathPbList->pixelBlockData, block.blockData);
                                                                }
                                                                // receivePixelBlockData(pathPbList->pixelBlockData, block.blockData);
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d wait the pixel %d Coord(%d,%d) for ev of process %d aluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }
                                                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                                {
                                                                    logW("process %d receive the pixel %d Coord(%d,%d) for eval of process %d uate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pathPbList->pixelBlockData.coreData.id, pathPbList->pixelBlockData.coreData.coord.i, pathPbList->pixelBlockData.coreData.coord.j, tabBlockData[pathPbList->pixelBlockData.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                                }
                                                            }
                                                        }
                                                    }
                                            }
                                        }
                                }
                                else
                                {
                                    pixel2 = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j - (pixel.coreData.id - id2P)];
                                }

                                if (id1P == id2P)
                                {
                                    if (pixel.coreData.id % 2 == 0)
                                    {
                                        begin = pixel1.coreData.firstBound.i - 1;
                                        end = begin;
                                    }
                                    else
                                    {
                                        begin = pixel1.coreData.firstBound.i - 1;
                                        end = pixel1.coreData.firstBound.j;
                                    }
                                    init = 1;
                                    if (d == 50000)
                                    {
                                        printPixelBlockData(pixel1);
                                        printPixelBlockData(pixel2);
                                        printf("begin %d end %d\n", begin, end);
                                    }
                                }
                                else
                                {
                                    begin = pixel1.coreData.firstBound.i - 1;
                                    end = pixel1.coreData.firstBound.j - 1;
                                    init = 0;
                                    if (d == 50000)
                                    {
                                        printPixelBlockData(pixel1);
                                        printPixelBlockData(pixel2);
                                        printf("begin %d end %d id1P = %d id2P = %d gap = %d\n", begin, end, id1P, id2P, gap);
                                    }
                                }

                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                {
                                    logD("process %d start the first (+,min) of the pixel %d Coord(%d,%d) (id1=%d; id2=%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, pixel1.coreData.id, pixel2.coreData.id, d, MPI_Wtime() - startTime);
                                }
                                tmpDouble = MPI_Wtime();
                                computeMCOP(pixel, begin, end, init, 2);
                                calculTime += MPI_Wtime() - tmpDouble;
                                if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                {
                                    logD("process %d finish the first (+,min) of the pixel %d Coord(%d,%d) (id1=%d; id2=%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, pixel1.coreData.id, pixel2.coreData.id, d, MPI_Wtime() - startTime);
                                }

                                if (id1P != id2P)
                                {
                                    if (id2P <= pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j]->pixelBlockData.coreData.id)
                                    {
                                        pathPbList = pixelatedBlock.needAllPixelBlockDataListByColumn.pixelBlockDataList[coord.j];
                                        while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id2P)
                                        {
                                            pathPbList = pathPbList->next;
                                        }
                                        if (pathPbList == NULL)
                                        {
                                            logE("process %d diagonal %d gap = %d => the list of need pixel of blocks must contains the id2P %d", rank, d, gap, id2P);
                                            exit(EXIT_FAILURE);
                                        }
                                        pixel1 = pathPbList->pixelBlockData;

                                        if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pixel1.idBlockData].rank)
                                            if (pixel1.idBlockData != block.blockData.coreData.id && (coord.i == 1)) //|| coord.j == 1
                                            {
                                                if (!isEmptyBlockDataList(block.needBlockDataList))
                                                {
                                                    pathBList = block.needBlockDataList;
                                                    while (pathBList != NULL && pixel1.idBlockData != pathBList->blockData.coreData.id)
                                                    {
                                                        pathBList = pathBList->next;
                                                    }
                                                    if (pathBList != NULL)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel1, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel1, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList4s(receivedList4s, pixel1.coreData.id, pixel1.coreData.address, pixel1.coreData.coord) == 0)
                                                        {
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d wait the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                            // checking before receiving
                                                            if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                            {
                                                                coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                {
                                                                    receivePixelBlockData4s(kBlock, block.blockData);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                receivePixelBlockData4s(pixel1, block.blockData);
                                                            }
                                                            // receivePixelBlockData(pixel1, block.blockData);
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d receive the pixel %d Coord(%d,%d) for  of process %d evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pixel1.coreData.id, pixel1.coreData.address, pixel1.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel1.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel1.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel1.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel1.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel1.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel1.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel1, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel1, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                    else
                                    {
                                        pixel1 = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i - (pixel.coreData.id - id2P)][coord.j];
                                    }

                                    if (id1P <= pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i]->pixelBlockData.coreData.id)
                                    {
                                        pathPbList = pixelatedBlock.needAllPixelBlockDataListByLine.pixelBlockDataList[coord.i];
                                        while (pathPbList != NULL && pathPbList->pixelBlockData.coreData.id != id1P)
                                        {
                                            pathPbList = pathPbList->next;
                                        }
                                        if (pathPbList == NULL)
                                        {
                                            logE("the list of need pixel of blocks must contains the id %d", id1P);
                                            exit(EXIT_FAILURE);
                                        }
                                        pixel2 = pathPbList->pixelBlockData;

                                        if (tabBlockData[pixel.idBlockData].rank != tabBlockData[pixel2.idBlockData].rank)
                                            if (pixel2.idBlockData != block.blockData.coreData.id && (coord.j == 1)) // coord.i == 1 ||
                                            {
                                                if (!isEmptyBlockDataList(block.needBlockDataList))
                                                {
                                                    pathBList = block.needBlockDataList;
                                                    while (pathBList != NULL && pixel2.idBlockData != pathBList->blockData.coreData.id)
                                                    {
                                                        pathBList = pathBList->next;
                                                    }
                                                    if (pathBList != NULL)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel2, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel2, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList4s(receivedList4s, pixel2.coreData.id, pixel2.coreData.address, pixel2.coreData.coord) == 0)
                                                        {
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d wait the pixel %d Coord(%d,%d) fo of process %d r evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                            // checking before receiving
                                                            if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                            {
                                                                coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                                PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                                if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                                {
                                                                    receivePixelBlockData4s(kBlock, block.blockData);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                receivePixelBlockData4s(pixel2, block.blockData);
                                                            }
                                                            // receivePixelBlockData(pixel2, block.blockData);
                                                            if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                            {
                                                                logW("process %d receive the pixel %d Coord(%d,%d) for  of process %d evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                            }
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList4s(receivedList4s, pixel2.coreData.id, pixel2.coreData.address, pixel2.coreData.coord) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        // checking before receiving
                                                        if (tabBlockData[pixel2.idBlockData].fragLevel < maxFrag - 1 && tabBlockData[pixel2.idBlockData].fragLevel != tabBlockData[pixel.idBlockData].fragLevel)
                                                        {
                                                            coordBis.i = (int)ceil(pixel2.coreData.coord.i / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            coordBis.j = (int)ceil(pixel2.coreData.coord.j / (float)pow(2, (tabBlockData[pixel.idBlockData].fragLevel - tabBlockData[pixel2.idBlockData].fragLevel - (tabBlockData[pixel.idBlockData].fragLevel == maxFrag ? 1 : 0))));
                                                            PixelBlockData kBlock = tabPixelatedBlockData[pixel2.idBlockData].pixelBlockDataTab[coordBis.i][coordBis.j];
                                                            if (findTagToReceivedList4s(receivedList4s, kBlock.coreData.id, kBlock.coreData.address, kBlock.coreData.coord) == 0)
                                                            {
                                                                receivePixelBlockData4s(kBlock, block.blockData);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            receivePixelBlockData4s(pixel2, block.blockData);
                                                        }
                                                        // receivePixelBlockData(pixel2, block.blockData);
                                                        if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d)  of process %d for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel2.coreData.id, pixel2.coreData.coord.i, pixel2.coreData.coord.j, tabBlockData[pixel2.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                    else
                                    {
                                        pixel2 = pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j - (pixel.coreData.id - id1P)];
                                    }

                                    begin = pixel2.coreData.secondBound.i;
                                    end = pixel2.coreData.secondBound.j;
                                    init = 0;
                                    if (d == 50000)
                                    {
                                        printPixelBlockData(pixel1);
                                        printPixelBlockData(pixel2);
                                        printf("begin %d end %d\n", begin, end);
                                    }
                                    /*printPixelBlockData(pixel1);
                                    printPixelBlockData(pixel2);
                                    printf("begin %d end %d\n", begin, end);*/

                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d start the second (+,min) of the pixel %d Coord(%d,%d) (id1=%d; id2=%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, pixel1.coreData.id, pixel2.coreData.id, d, MPI_Wtime() - startTime);
                                    }
                                    tmpDouble = MPI_Wtime();
                                    computeMCOP(pixel, begin, end, init, 2);
                                    calculTime += MPI_Wtime() - tmpDouble;
                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d finish the second (+,min) of the pixel %d Coord(%d,%d) (id1=%d; id2=%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, pixel1.coreData.id, pixel2.coreData.id, d, MPI_Wtime() - startTime);
                                    }
                                }

                                if (finalize)
                                {

                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d start the finalisation of the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                    }

                                    tmpDouble = MPI_Wtime();
                                    computeMCOP(pixel, 0, 0, 0, 3);
                                    calculTime += MPI_Wtime() - tmpDouble;

                                    if (d == 3 && verboseDebug && rank == 7) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d finish the finalisation of the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                    }

                                    tmpDouble = MPI_Wtime();
                                    BlockDataList *list1 = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                                    list1 = deleteOccurBlockList(list1);
                                    BlockDataList *list2 = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                                    list2 = deleteOccurBlockList(list2);

                                    BlockDataList *path = buildSendBlockDataList(list1, list2);

                                    if (!isEmptyBlockDataList(path))
                                        if (d == 10000 && verboseDebug) // && (pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 2 || pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id == 3))
                                        {
                                            logE("process %d send the pixel %d Coord(%d,%d) d = %d to process %d in %f s \t==> a = %d nbEvaluatePixel = %d diag = %d", rank, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j].coreData.id, coord.i, coord.j, d, path->blockData.rank, MPI_Wtime() - startTime, a, pixelatedBlock.pBlockData.nbEvaluatePixel, pixelatedBlock.pBlockData.maxDiag);
                                        }
                                    if (a == pixelatedBlock.pBlockData.maxDiag)
                                        sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 1);
                                    else
                                        sendPixelBlockData(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path, 0);
                                    buildTime += MPI_Wtime() - tmpDouble;
                                }
                            }
                        }

                        if (a % 2 != 0)
                            gap--;

                        if (finalize == 1)
                        {
                            if (a == pixelatedBlock.pBlockData.maxDiag)
                                stop = 1;
                            else
                            {
                                step++;
                                if (a % 2 != 0)
                                    tmpGap--;
                                finalize = 0;
                            }
                        }
                    }

                    if (stop != 1)
                    {
                        if (diagEval != pixelatedBlock.pBlockData.maxDiag)
                            diagEval += 2;
                        gap = tmpGap;
                        tmpGap = ++gap;
                    }
                }
                // printf("Finish diag %ld => rank %d\n", d, rank);
            }

            if ((d == 1 || d == 2 || d == 3 || d == 4 || d == 5) && block.blockData.rank == rank && block.blockData.coreData.diag == d && verboseDebug)
            {
                logI("process %d finishs d = %d in %f s", rank, d, MPI_Wtime() - startTime);
            }

            i++;
            if (i != maxEvalBlock)
            {
                block = tabBlock[i];
                pixelatedBlock = tabPixelatedBlock[i];
            }
        }

        if ((rank == 4 || rank == 5 || rank == 6) && d == 200)
        {
            logD("process %d finishs d = %d in %f s", rank, d, MPI_Wtime() - startTime);
        }

        if ((rank == 5 || rank == 6 || rank == 7) && d == 500)
        {
            logD("process %d finishs d = %d in %f s", rank, d, MPI_Wtime() - startTime);
        }

        if ((rank == 0 || rank == 3) && d == 200)
        {
            logD("process %d finishs d = %d in %f s", rank, d, MPI_Wtime() - startTime);
        }

        if ((rank == 1 || rank == 2 || rank == 3) && d == 400)
        {
            logD("process %d finishs d = %d in %f s", rank, d, MPI_Wtime() - startTime);
        }

        /*if ((d == 2 && (rank == 6 || rank == 7 || rank == 8 || rank == 9 || rank == 10)) || ((d == 3 && (rank == 11 || rank == 12 || rank == 13 || rank == 14))))
        {
            logD("process %d finishs d = %d in %f s", rank, d, MPI_Wtime() - startTime);
        }*/
    }
    // free(receivedList4s);

    if (last == 1000 || rank == 40000)
    {
        int i = 0;
        int j = 0;
        FILE *file = fopen("test_p0.csv", "a");
        printf("**********************************************************************\n");
        printf("*                       Dynamic programming Table                    *\n");
        printf("**********************************************************************\n");
        for (i = 1; i <= maxNumber; i++)
        {
            for (j = 1; j <= maxNumber; j++)
            {
                if (i > j)
                    // printf("\t");
                    fprintf(file, ";");
                else
                    // printf("%d\t", getMCOP(i, j));
                    fprintf(file, "%d;", getMCOP(i, j));
            }
            // printf("\n");
            fprintf(file, "\n");
        }
    }

    // logD("rank %d : build time ==> %f s\n", rank, buildTime);
    return (last ? getMCOP(1, maxNumber) : 0);
}

void receivePixelBlockData4s(PixelBlockData src, BlockData dest)
{
    int count, blockLengths, stride, absc, ord;

    count = src.coreData.dim.nbLine;
    blockLengths = src.coreData.dim.nbColumn;
    stride = maxNumber + 1;
    MPI_Type_vector(count, blockLengths, stride, MPI_INT, &vector);
    MPI_Type_commit(&vector);
    absc = src.coreData.firstBound.i;
    ord = src.coreData.secondBound.i;

    int flag;
    int tag = src.coreData.address;
    int source = tabBlockData[src.idBlockData].rank;

    MPI_Iprobe(source, tag, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);

    if (flag)
    {
        tmpDouble = MPI_Wtime();
        MPI_Recv(&tabMCOP[absc][ord], 1, vector, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        comTime += MPI_Wtime() - tmpDouble;
    }
    else
    {
        // MPI_Request request;
        tmpDouble = MPI_Wtime();
        MPI_Send(&tag, 1, MPI_INT, source, tag, MPI_COMM_WORLD); //, &request);

        MPI_Recv(&tabMCOP[absc][ord], 1, vector, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        comTime += MPI_Wtime() - tmpDouble;
    }

    if (rank == 10000 && tag == 960000)
    {
        printf("COMPIIIIIIIIIIIIIII\n");
        printPixelBlockData(src);
    }
    receivedList4s = addTagToReceivedList4s(receivedList4s, src.coreData.id, tag, src.coreData.coord);
    MPI_Type_free(&vector);
}

ReceivedList4s *createReceivedList4s()
{
    ReceivedList4s *receivedList4s = (ReceivedList4s *)malloc(sizeof *receivedList4s);
    receivedList4s->id = -1;
    receivedList4s->tag = -1;
    receivedList4s->coord.i = -1;
    receivedList4s->coord.j = -1;
    receivedList4s->next = NULL;
    return receivedList4s;
}
ReceivedList4s *addTagToReceivedList4s(ReceivedList4s *receivedList4s, int id, int tag, Coord coord)
{
    if (isEmptyReceivedList4s(receivedList4s))
    {
        receivedList4s->id = id;
        receivedList4s->tag = tag;
        receivedList4s->coord = coord;
        return receivedList4s;
    }
    else
    {
        ReceivedList4s *new = createReceivedList4s();
        new->id = id;
        new->tag = tag;
        new->coord = coord;
        new->next = receivedList4s;
        return new;
    }
}
int isEmptyReceivedList4s(ReceivedList4s *receivedList4s) { return receivedList4s->tag == -1 && receivedList4s->id == -1; }
int findTagToReceivedList4s(ReceivedList4s *receivedList4s, int id, int tag, Coord coord)
{
    if (isEmptyReceivedList4s(receivedList4s))
        return 0;
    ReceivedList4s *path = receivedList4s;
    while (path != NULL && !(path->id == id && path->tag == tag && path->coord.i == coord.i && path->coord.j == coord.j))
        path = path->next;
    return (path != NULL ? 1 : 0);
}