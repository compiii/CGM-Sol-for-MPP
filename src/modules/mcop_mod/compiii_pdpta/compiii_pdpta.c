/**
 * @file compiii_pdpta.c
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief 
 * @version 1.2
 * @date 2020-07-17
 * 
 * @copyright Compiii Thesis Copyright (c) 2020
 * 
 */

#include "compiii_pdpta.h"
#include "clogger.h"

int compiiiPdpta()
{
    long d = 0;
    int i = 0, last = 0, a, b, iter;
    Block block;
    PixelatedBlock pixelatedBlock;
    Coord coord;
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
                            tmpDouble = MPI_Wtime();
                            computeMCOP(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], 0, 0, 0, 0);
                            calculTime += MPI_Wtime() - tmpDouble;
                        }
                    }
                    else
                    {
                        iter = 0;
                        for (b = 2 + (a - pixelatedBlock.pBlockData.nbEvaluatePixel - 1); b <= pixelatedBlock.pBlockData.nbEvaluatePixel; b++)
                        {
                            coord.i = pixelatedBlock.pBlockData.nbEvaluatePixel - iter++;
                            coord.j = b;
                            tmpDouble = MPI_Wtime();
                            computeMCOP(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], 0, 0, 0, 1);
                            calculTime += MPI_Wtime() - tmpDouble;
                        }
                    }
                }
                for (a = pixelatedBlock.pBlockData.nbEvaluatePixel; a <= pixelatedBlock.pBlockData.maxDiag; a++)
                {
                    if (a <= pixelatedBlock.pBlockData.nbEvaluatePixel)
                    {
                        for (b = 1; b <= a; b++)
                        {
                            coord.i = a - b + 1;
                            coord.j = b;
                            BlockDataList *path = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                            sendedList = createReceivedList();
                            if (!isEmptyBlockDataList(path))
                                while (path != NULL)
                                {
                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                    }
                                    path = path->next;
                                }

                            free(sendedList);
                            sendedList = createReceivedList();
                            path = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                            if (!isEmptyBlockDataList(path))
                                while (path != NULL)
                                {
                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                    }
                                    path = path->next;
                                }
                            free(sendedList);
                        }
                    }
                    else
                    {
                        iter = 0;
                        for (b = 2 + (a - pixelatedBlock.pBlockData.nbEvaluatePixel - 1); b <= pixelatedBlock.pBlockData.nbEvaluatePixel; b++)
                        {
                            coord.i = pixelatedBlock.pBlockData.nbEvaluatePixel - iter++;
                            coord.j = b;
                            BlockDataList *path = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                            sendedList = createReceivedList();
                            if (!isEmptyBlockDataList(path))
                                while (path != NULL)
                                {
                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                    }
                                    path = path->next;
                                }

                            free(sendedList);
                            sendedList = createReceivedList();
                            path = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                            if (!isEmptyBlockDataList(path))
                                while (path != NULL)
                                {
                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                    }
                                    path = path->next;
                                }
                            free(sendedList);
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
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(0, path->pixelBlockData, block.blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        path = path->next;
                                    }
                            }

                            if (coord.j == 1)
                            {
                                PixelBlockDataList *path = pixelatedBlock.needPixelBlockDataListByLine.pixelBlockDataList[coord.i];
                                if (!isEmptyPixelBlockDataList(path))
                                    while (path != NULL)
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(0, path->pixelBlockData, block.blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        path = path->next;
                                    }
                            }

                            tmpDouble = MPI_Wtime();
                            computeMCOP(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], 0, 0, 0, 0);
                            calculTime += MPI_Wtime() - tmpDouble;
                        }
                    }
                    else
                    {
                        iter = 0;
                        for (b = 2 + (a - pixelatedBlock.pBlockData.nbEvaluatePixel - 1); b <= pixelatedBlock.pBlockData.nbEvaluatePixel; b++)
                        {
                            coord.i = pixelatedBlock.pBlockData.nbEvaluatePixel - iter++;
                            coord.j = b;

                            tmpDouble = MPI_Wtime();
                            computeMCOP(pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], 0, 0, 0, 1);
                            calculTime += MPI_Wtime() - tmpDouble;
                        }
                    }
                }
                for (a = 1; a <= pixelatedBlock.pBlockData.maxDiag; a++)
                {
                    if (a <= pixelatedBlock.pBlockData.nbEvaluatePixel)
                    {
                        for (b = 1; b <= a; b++)
                        {
                            coord.i = a - b + 1;
                            coord.j = b;
                            BlockDataList *path = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                            sendedList = createReceivedList();
                            if (!isEmptyBlockDataList(path))
                                while (path != NULL)
                                {
                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                    }
                                    path = path->next;
                                }

                            free(sendedList);
                            sendedList = createReceivedList();
                            path = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                            if (!isEmptyBlockDataList(path))
                                while (path != NULL)
                                {
                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                    }
                                    path = path->next;
                                }
                            free(sendedList);
                        }
                    }
                    else
                    {
                        iter = 0;
                        for (b = 2 + (a - pixelatedBlock.pBlockData.nbEvaluatePixel - 1); b <= pixelatedBlock.pBlockData.nbEvaluatePixel; b++)
                        {
                            coord.i = pixelatedBlock.pBlockData.nbEvaluatePixel - iter++;
                            coord.j = b;
                            BlockDataList *path = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                            sendedList = createReceivedList();
                            if (!isEmptyBlockDataList(path))
                                while (path != NULL)
                                {
                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                    }
                                    path = path->next;
                                }

                            free(sendedList);
                            sendedList = createReceivedList();
                            path = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                            if (!isEmptyBlockDataList(path))
                                while (path != NULL)
                                {
                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                    {
                                        tmpDouble = MPI_Wtime();
                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                        comTime += MPI_Wtime() - tmpDouble;
                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                    }
                                    path = path->next;
                                }
                            free(sendedList);
                        }
                    }
                }
            }
            else
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

                                if (rank == 3000 && d == 3)
                                    logD("id1P = %d id2P = %d gap = %d\n", id1P, id2P, gap);

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
                                                    if (rank == 1100 && d == 3)
                                                    {
                                                        logD("id1P = %d id2P = %d gap = %d\n", id1P, id2P, gap);
                                                        printPixelBlockData(pixel1);
                                                    }

                                                    if (d == 3 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                    tmpDouble = MPI_Wtime();
                                                    communicateBlockData(0, pixel1, block.blockData);
                                                    comTime += MPI_Wtime() - tmpDouble;
                                                    if (d == 3 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pixel1.coreData.address) == 0)
                                                    {
                                                        if (d == 3 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel1, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                        if (d == 3 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                        {
                                                            logW("process %d receive the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (findTagToReceivedList(receivedList, pixel1.coreData.address) == 0)
                                                {
                                                    if (d == 3 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d wait the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                                    }
                                                    tmpDouble = MPI_Wtime();
                                                    communicateBlockData(0, pixel1, block.blockData);
                                                    comTime += MPI_Wtime() - tmpDouble;
                                                    if (d == 3 && verboseDebug) // && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                                    {
                                                        logW("process %d receive the pixel %d Coord(%d,%d) of process %d  for evaluate the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel1.coreData.id, pixel1.coreData.coord.i, pixel1.coreData.coord.j, tabBlockData[pixel1.idBlockData].rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
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
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList(receivedList, pathPbList->pixelBlockData.coreData.address) == 0)
                                                        {
                                                            tmpDouble = MPI_Wtime();
                                                            communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                            comTime += MPI_Wtime() - tmpDouble;
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pathPbList->pixelBlockData.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
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
                                        if (pixel2.idBlockData != block.blockData.coreData.id && (coord.j == 1)) //coord.i == 1 ||
                                        {
                                            if (!isEmptyBlockDataList(block.needBlockDataList))
                                            {
                                                pathBList = block.needBlockDataList;
                                                while (pathBList != NULL && pixel2.idBlockData != pathBList->blockData.coreData.id)
                                                {
                                                    if (rank == 200 && d == 4)
                                                    {
                                                        logD("id1P = %d id2P = %d pixel = %d pixel2 = %d gap = %d\n", id1P, id2P, pixel.coreData.id, pixel2.coreData.id, gap);
                                                        //printPixelBlockData(pathPbList->pixelBlockData);
                                                    }
                                                    if (rank == 200 && d == 4)
                                                        printBlockData(pathBList->blockData);
                                                    pathBList = pathBList->next;
                                                }
                                                if (pathBList != NULL)
                                                {
                                                    tmpDouble = MPI_Wtime();
                                                    communicateBlockData(0, pixel2, block.blockData);
                                                    comTime += MPI_Wtime() - tmpDouble;
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pixel2.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel2, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (findTagToReceivedList(receivedList, pixel2.coreData.address) == 0)
                                                {
                                                    tmpDouble = MPI_Wtime();
                                                    communicateBlockData(0, pixel2, block.blockData);
                                                    comTime += MPI_Wtime() - tmpDouble;
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
                                            if (pathPbList->pixelBlockData.idBlockData != block.blockData.coreData.id && (coord.j == 1)) //coord.i == 1 ||
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
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList(receivedList, pathPbList->pixelBlockData.coreData.address) == 0)
                                                        {
                                                            tmpDouble = MPI_Wtime();
                                                            communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                            comTime += MPI_Wtime() - tmpDouble;
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pathPbList->pixelBlockData.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
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
                                    if (d == 50000)
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

                                if (d == 3 && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                {
                                    logD("process %d start the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                }

                                tmpDouble = MPI_Wtime();
                                computeMCOP(pixel, begin, end, init, 2);
                                calculTime += MPI_Wtime() - tmpDouble;

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
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel1, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList(receivedList, pixel1.coreData.address) == 0)
                                                        {
                                                            tmpDouble = MPI_Wtime();
                                                            communicateBlockData(0, pixel1, block.blockData);
                                                            comTime += MPI_Wtime() - tmpDouble;
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pixel1.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel1, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
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
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel2, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList(receivedList, pixel2.coreData.address) == 0)
                                                        {
                                                            tmpDouble = MPI_Wtime();
                                                            communicateBlockData(0, pixel2, block.blockData);
                                                            comTime += MPI_Wtime() - tmpDouble;
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pixel2.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel2, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
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

                                    tmpDouble = MPI_Wtime();
                                    computeMCOP(pixel, begin, end, init, 2);
                                    calculTime += MPI_Wtime() - tmpDouble;
                                }

                                if (finalize)
                                {
                                    tmpDouble = MPI_Wtime();
                                    computeMCOP(pixel, 0, 0, 0, 3);
                                    calculTime += MPI_Wtime() - tmpDouble;

                                    if (d == 3 && (pixel.coreData.id == 4 || pixel.coreData.id == 6))
                                    {
                                        logD("process %d finishs the pixel %d Coord(%d,%d) d = %d in %f s", rank, pixel.coreData.id, pixel.coreData.coord.i, pixel.coreData.coord.j, d, MPI_Wtime() - startTime);
                                    }
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
                                                    tmpDouble = MPI_Wtime();
                                                    communicateBlockData(0, pixel1, block.blockData);
                                                    comTime += MPI_Wtime() - tmpDouble;
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pixel1.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel1, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (findTagToReceivedList(receivedList, pixel1.coreData.address) == 0)
                                                {
                                                    tmpDouble = MPI_Wtime();
                                                    communicateBlockData(0, pixel1, block.blockData);
                                                    comTime += MPI_Wtime() - tmpDouble;
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
                                                                tmpDouble = MPI_Wtime();
                                                                communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                                comTime += MPI_Wtime() - tmpDouble;
                                                            }
                                                            else
                                                            {
                                                                if (findTagToReceivedList(receivedList, pathPbList->pixelBlockData.coreData.address) == 0)
                                                                {
                                                                    tmpDouble = MPI_Wtime();
                                                                    communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                                    comTime += MPI_Wtime() - tmpDouble;
                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                            if (findTagToReceivedList(receivedList, pathPbList->pixelBlockData.coreData.address) == 0)
                                                            {
                                                                tmpDouble = MPI_Wtime();
                                                                communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                                comTime += MPI_Wtime() - tmpDouble;
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
                                        if (pixel2.idBlockData != block.blockData.coreData.id && (coord.j == 1)) //coord.i == 1 ||
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
                                                    tmpDouble = MPI_Wtime();
                                                    communicateBlockData(0, pixel2, block.blockData);
                                                    comTime += MPI_Wtime() - tmpDouble;
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pixel2.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel2, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (findTagToReceivedList(receivedList, pixel2.coreData.address) == 0)
                                                {
                                                    tmpDouble = MPI_Wtime();
                                                    communicateBlockData(0, pixel2, block.blockData);
                                                    comTime += MPI_Wtime() - tmpDouble;
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
                                                    if (pathPbList->pixelBlockData.idBlockData != block.blockData.coreData.id && (coord.j == 1)) //coord.i == 1 ||
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
                                                                tmpDouble = MPI_Wtime();
                                                                communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                                comTime += MPI_Wtime() - tmpDouble;
                                                            }
                                                            else
                                                            {
                                                                if (findTagToReceivedList(receivedList, pathPbList->pixelBlockData.coreData.address) == 0)
                                                                {
                                                                    tmpDouble = MPI_Wtime();
                                                                    communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                                    comTime += MPI_Wtime() - tmpDouble;
                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                            if (findTagToReceivedList(receivedList, pathPbList->pixelBlockData.coreData.address) == 0)
                                                            {
                                                                tmpDouble = MPI_Wtime();
                                                                communicateBlockData(0, pathPbList->pixelBlockData, block.blockData);
                                                                comTime += MPI_Wtime() - tmpDouble;
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

                                tmpDouble = MPI_Wtime();
                                computeMCOP(pixel, begin, end, init, 2);
                                calculTime += MPI_Wtime() - tmpDouble;

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
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel1, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList(receivedList, pixel1.coreData.address) == 0)
                                                        {
                                                            tmpDouble = MPI_Wtime();
                                                            communicateBlockData(0, pixel1, block.blockData);
                                                            comTime += MPI_Wtime() - tmpDouble;
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pixel1.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel1, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
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
                                            if (pixel2.idBlockData != block.blockData.coreData.id && (coord.j == 1)) //coord.i == 1 ||
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
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel2, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                    }
                                                    else
                                                    {
                                                        if (findTagToReceivedList(receivedList, pixel2.coreData.address) == 0)
                                                        {
                                                            tmpDouble = MPI_Wtime();
                                                            communicateBlockData(0, pixel2, block.blockData);
                                                            comTime += MPI_Wtime() - tmpDouble;
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    if (findTagToReceivedList(receivedList, pixel2.coreData.address) == 0)
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(0, pixel2, block.blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
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

                                    tmpDouble = MPI_Wtime();
                                    computeMCOP(pixel, begin, end, init, 2);
                                    calculTime += MPI_Wtime() - tmpDouble;
                                }

                                if (finalize)
                                {
                                    tmpDouble = MPI_Wtime();
                                    computeMCOP(pixel, 0, 0, 0, 3);
                                    calculTime += MPI_Wtime() - tmpDouble;
                                }
                            }
                        }

                        if (a % 2 != 0)
                            gap--;

                        if (finalize == 1)
                        {
                            if (a == pixelatedBlock.pBlockData.maxDiag)
                            {
                                stop = 1;
                                for (a = 1; a <= pixelatedBlock.pBlockData.maxDiag; a++)
                                {
                                    if (a <= pixelatedBlock.pBlockData.nbEvaluatePixel)
                                    {
                                        for (b = 1; b <= a; b++)
                                        {
                                            coord.i = a - b + 1;
                                            coord.j = b;
                                            BlockDataList *path = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                                            sendedList = createReceivedList();
                                            if (!isEmptyBlockDataList(path))
                                                while (path != NULL)
                                                {
                                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                                    }
                                                    path = path->next;
                                                }

                                            free(sendedList);
                                            sendedList = createReceivedList();
                                            path = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                                            if (!isEmptyBlockDataList(path))
                                                while (path != NULL)
                                                {
                                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                                    }
                                                    path = path->next;
                                                }
                                            free(sendedList);
                                        }
                                    }
                                    else
                                    {
                                        iter = 0;
                                        for (b = 2 + (a - pixelatedBlock.pBlockData.nbEvaluatePixel - 1); b <= pixelatedBlock.pBlockData.nbEvaluatePixel; b++)
                                        {
                                            coord.i = pixelatedBlock.pBlockData.nbEvaluatePixel - iter++;
                                            coord.j = b;
                                            BlockDataList *path = pixelatedBlock.dependBlockDataListByLine.blockDataList[coord.i];
                                            sendedList = createReceivedList();
                                            if (!isEmptyBlockDataList(path))
                                                while (path != NULL)
                                                {
                                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                                    }
                                                    path = path->next;
                                                }

                                            free(sendedList);
                                            sendedList = createReceivedList();
                                            path = pixelatedBlock.dependBlockDataListByColumn.blockDataList[coord.j];
                                            if (!isEmptyBlockDataList(path))
                                                while (path != NULL)
                                                {
                                                    if (!findTagToReceivedList(sendedList, path->blockData.rank))
                                                    {
                                                        tmpDouble = MPI_Wtime();
                                                        communicateBlockData(1, pixelatedBlock.pBlockData.pixelBlockDataTab[coord.i][coord.j], path->blockData);
                                                        comTime += MPI_Wtime() - tmpDouble;
                                                        sendedList = addTagToReceivedList(sendedList, path->blockData.rank);
                                                    }
                                                    path = path->next;
                                                }
                                            free(sendedList);
                                        }
                                    }
                                }
                            }
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
                //printf("Finish diag %ld => rank %d\n", d, rank);
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

        if ((d == 2 && (rank == 6 || rank == 7 || rank == 8 || rank == 9 || rank == 10)) || ((d == 3 && (rank == 11 || rank == 12 || rank == 13 || rank == 14))))
        {
            logD("process %d finishs d = %d in %f s", rank, d, MPI_Wtime() - startTime);
        }
    }
    free(receivedList);

    if (last == 5000 || rank == 200000)
    {
        int i = 0;
        int j = 0;
        printf("**********************************************************************\n");
        printf("*                       Dynamic programming Table                    *\n");
        printf("**********************************************************************\n");
        for (i = 1; i <= maxNumber; i++)
        {
            for (j = 1; j <= maxNumber; j++)
            {
                if (i > j)
                    printf("\t");
                else
                    printf("%d\t", getMCOP(i, j));
            }
            printf("\n");
        }
    }

    return (last ? getMCOP(1, maxNumber) : 0);
}