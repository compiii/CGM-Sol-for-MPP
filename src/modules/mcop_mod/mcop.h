/**
 * @file mcop.h
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief 
 * @version 1.1
 * @date 2019-11-05
 * 
 * @copyright Compiii Thesis Copyright (c) 2019
 * 
 */
#ifndef MCOP_H_
#define MCOP_H_

#include "main.h"

typedef struct SentList_T SentList;
struct SentList_T
{
    PixelBlockData src;
    BlockData dest;
    SentList *next;
};

SentList *sentList;

SentList *createSentList();
int isEmptySentList(SentList *list);
SentList *addDataToSentList(SentList *list, PixelBlockData src, BlockData dest);
BlockDataList *removeBlockDataToBlockDataList(BlockDataList *list, BlockData blockData);
BlockDataList *deleteOccurBlockList(BlockDataList *list);
BlockDataList *buildSendBlockDataList(BlockDataList *list1, BlockDataList *list2);

void sendPixelBlockData(PixelBlockData src, BlockDataList *dependBlockList, int finalized);
void receivePixelBlockData(PixelBlockData src, BlockData dest);

typedef struct ReceivedList_T ReceivedList;
struct ReceivedList_T
{
    int tag;
    ReceivedList *next;
};

ReceivedList *receivedList, *sendedList;

int **tabMCOP;

void initMCOPMod();
int godboleAlgorithm();
int getMCOP(int i, int j);
int max(int n, int m);
int min(int m, int n);

void computeMCOP(PixelBlockData pixel, int begin, int end, int init, int option);
void communicateBlockData(int type, PixelBlockData src, BlockData blockData);
void receiveBlockData(BlockData blockData, MPI_Datatype vector, unsigned long absc, unsigned long ord, int tag);
void sendBlockData(BlockData blockData, MPI_Datatype vector, unsigned long absc, unsigned long ord, int tag);

ReceivedList *createReceivedList();
ReceivedList *addTagToReceivedList(ReceivedList *receivedList, int tag);
int isEmptyReceivedList(ReceivedList *receivedList);
int findTagToReceivedList(ReceivedList *receivedList, int tag);

#endif
