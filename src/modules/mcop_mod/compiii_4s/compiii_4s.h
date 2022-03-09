/**
 * @file compiii_4s.h
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-01-10
 *
 * @copyright Compiii Thesis Copyright (c) 2022
 *
 */
#ifndef COMPIII_4S
#define COMPIII_4S

#include "mcop.h"

typedef struct ReceivedList4s_T ReceivedList4s;
struct ReceivedList4s_T
{
    int id;
    int tag;
    Coord coord;
    ReceivedList4s *next;
};

ReceivedList4s *receivedList4s;

int compiii4s();
void receivePixelBlockData4s(PixelBlockData src, BlockData dest);

ReceivedList4s *createReceivedList4s();
ReceivedList4s *addTagToReceivedList4s(ReceivedList4s *receivedList, int id, int tag, Coord coord);
int isEmptyReceivedList4s(ReceivedList4s *receivedList);
int findTagToReceivedList4s(ReceivedList4s *receivedList, int id, int tag, Coord coord);

#endif // !COMPIII_4S