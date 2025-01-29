//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// $Id$

#include <limits>
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <fcntl.h>

#include "swf.h"

#define SWFVERSION 3

using namespace std;

/** \fn void bi(char* a, int nbytes)
 * \brief Convert a string of char from little-endian to big-endian (byte order)
 * \param[in,out] a Character string that starts out as little-endian
 * \param[in] nbytes number of bytes to convert
 */
void
bi(char* a, int nbytes);

/** \fn int16_t nsb(int16_t x)
 * \brief Returns the min number of bits necessary to hold signed 16-bit
 * int "x"
 * \param[in] x A 16 bit integer.
 */
int16_t
nsb(int16_t x);

/** \fn int bytealign(int nbits)
 * \brief Returns the number of bytes required to hold nbits bits
 */
int
bytealign(int nbits);


/** \fn void copyBits(int nb, char* b, int bstart, char* a, int astart)
 * \brief Coby "nb" bits from "b" starting at "bstart" to "a" starting at bit
 * "astart"
 * \param[in] nb number of bits to copy from b
 * \param[in] b Array to copy from.
 * \param[in] bstart bit index in b to start copying from.
 * \param[in,out] a buffer to copy bits into 
 * \param[in] astart bit index in a to start copying to
 */
void
copyBits (int nb, char* b, int bstart, char* a, int astart);

/////////////////////////////  Helper Functions ////////////////////////////////
void
bi(char* a, int nbytes) {
#ifndef WORDS_BIGENDIAN
	int j = nbytes-1;
	for (int i=0; i<nbytes/2; i++) {
		char c = a[i];
		a[i] = a[j];
		a[j] = c;
		j--;
	}
#endif
}

int16_t
nsb(int16_t x) {
// Returns the min number of bits necessary to
// hold signed 16-bit int "x"
	if (x < 0)
		x = -x;
	double log2x = log10((double)x)/log10(2.0);
	int16_t rval = (int16_t)log2x;
	// two extra: one for the sign, one for good luck
	rval += 2;
	return rval;
}

int
bytealign(int nbits) {
// How many bytes does it take to hold "nbits" bits?
	int rval = nbits/8;
	if (nbits%8)
		rval++;
	return rval;
}

void
copyBits (int nb, char* b, int bstart, char* a, int astart) {
/*
 * Copy "nb" bits from "b" starting with "bstart"
 * to "a" starting at bit "astart"
 */
	int ia = astart;
	int ib = bstart;
	for (int i=0; i<nb; i++) {
		int ja = ia/8;		// byte number in a
		int jb = ib/8;		// byte number in b
		int ka = 7 - ia%8;  // bit number from low-bit
		int kb = 7 - ib%8;  // bit number from low-bit
		a[ja] |= ((((b[jb] >> kb)) & 1) << ka);
		ia++;
		ib++;
	}
}

////////////////////////// swf class ///////////////////////////////////////////


void 
swf::
makeRecordHeader(int code, int length) {
	RecordHeader temp(code, length);
	myNBytes += temp.write(&myBuf[myNBytes]);
}

void 
swf::
makeBackgroundColor( int color[3]) {
	makeRecordHeader(9, 3);

	myBuf[myNBytes++] = color[0];
	myBuf[myNBytes++] = color[1];
	myBuf[myNBytes++] = color[2];
}

void
swf::
makeFileHeader(int xmin, int xmax, int ymin, int ymax,
			   uint16_t nframes) {

	myBuf[myNBytes++] = 'F';
	myBuf[myNBytes++] = 'W';
	myBuf[myNBytes++] = 'S';
	myBuf[myNBytes++] = SWFVERSION;

	// uint32_t placeholder: write actual file size when done
	myBuf[myNBytes++] = 0;
	myBuf[myNBytes++] = 0;
	myBuf[myNBytes++] = 0;
	myBuf[myNBytes++] = 0;

	// frame size: xmin, xmax, ymin, ymax are in twips
	makeRect(xmin, xmax, ymin, ymax);

	// frame delay: high-order byte is left of decimal, low-order
	// is to the right (but is ignored - see pg 257) 12 fps
	myBuf[myNBytes++] = 0;	// low-order byte
	myBuf[myNBytes++] = (char)myFps;	// high-order byte

	// frame count
	memcpy(&myBuf[myNBytes], (char*)&nframes, sizeof(uint16_t)); // little-endian ui16
	myNBytes += sizeof(uint16_t);
}

void 
swf::
makeRect(int xmin, int xmax, int ymin, int ymax) {

	Rect temp(xmin, xmax, ymin, ymax);
	myNBytes += temp.write(&myBuf[myNBytes]);
}



//change from opengl coords to shockwave coords
void 
swf::
changeCoords(vector<FlashPoint> & points) {

	//scale and shift
	for (size_t i=0; i<points.size(); i++) {
		points[i].xC = myScale*(points[i].xC - myXShift);
		points[i].yC = myScale*(points[i].yC - myYShift);
	}

	// invert
	for (size_t i=0; i<points.size(); i++) {
		points[i].yC = myYMax - points[i].yC;
	}
}

//overloaded function for lines
void 
swf::
changeCoords(vector<FlashLine> & lines) {

	//scale and shift
	for (size_t i=0; i<lines.size(); i++) {
		lines[i].x[0] = myScale*(lines[i].x[0] - myXShift);
		lines[i].y[0] = myScale*(lines[i].y[0] - myYShift);

		lines[i].x[1] = myScale*(lines[i].x[1] - myXShift);
		lines[i].y[1] = myScale*(lines[i].y[1] - myYShift);
	}
	//invert
	for (size_t i=0; i<lines.size(); i++) {
		lines[i].y[0] = myYMax - lines[i].y[0];
		lines[i].y[1] = myYMax - lines[i].y[1];
	}
}

swf::
swf() : 
myBuf(1000000,0)
{
	myNBytes = 0;

	myId = 1;
	myTime = 4;
	myNframes = 40;
	myFps = myNframes/myTime;
	myXShift = 0;
	myYShift = 0;
	myScale = 1;
	myYMax = 0;

	//	myNumLineBits = 2;
	//	myNumFillBits = 1;

	myFramePixels = 300;

	myBGColor[0] = 255;
	myBGColor[1] = 255;
	myBGColor[2] = 255;

}

void 
swf::
setScale(int scaleVal){
	myScale = scaleVal;
}

void
swf::
setBackgroundColor(int color[]){
	myBGColor[0] = color[0];
	myBGColor[1] = color[1];
	myBGColor[2] = color[2];
}

void 
swf::
setNframes(int numFrames) {
	myNframes = numFrames;
	myFps = myNframes / myTime;
}

void 
swf::
setFramePixels(int pixelCount) {
	myFramePixels = pixelCount;
}

void 
swf::
setFrameDimen(int frameWidth, int frameHeight) {
	myFrameWidth = frameWidth;
	myFrameHeight = frameHeight;
}

void
swf::
setTime(int seconds) {
	myTime = seconds;
	myFps = myNframes / myTime;
}

void 
swf::
setShift(int XShift, int YShift) {
	myXShift = XShift;
	myYShift = YShift;
}

void
swf::
setYMax(int YMax) {
	myYMax = YMax;
}

void
swf::
openBuffer() {
	// determine framewidth and frameheight based on 
	// framepixels and scale
	makeFileHeader(0, myFrameWidth, 0, myFrameHeight, myNframes);

#if SWFVERSION > 7
	makeRecordHeader(69, 32);
	buf[nbytes++] = 0;
	buf[nbytes++] = 0;
	buf[nbytes++] = 0;
	buf[nbytes++] = 0;
#endif

	makeBackgroundColor(myBGColor);
}

void
swf::
closeBuffer() {
	// End tag - no body, just tag = 0
	makeRecordHeader(0, 0);
	memcpy(&myBuf[4], &myNBytes, 4);
}

void
swf::
writeBuffer(string path){
	mode_t perm = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
	int fd = open(path.c_str(), O_WRONLY | O_CREAT, perm);
	write(fd, &myBuf[0], myNBytes);
	close(fd);
}

void 
swf::
makeShape( const vector<GLfloat> & buffer) {

	vector<FlashLine> lines;
	vector<FlashPoint> points;

	parseFeedBuffer(buffer, lines, points);
	changeCoords(lines);
	changeCoords(points);
	DefineShape aShape(myId, lines, points);
	if (myNBytes + aShape.size() >= myBuf.size() ) {
	    myBuf.resize(myBuf.size()*2,0);
	}

	myNBytes += aShape.write(&myBuf[myNBytes]);
	myId++;

}

void
swf::
removeObject(int depth) {
	makeRecordHeader(28, 2);
	memcpy (&myBuf[myNBytes], (char*)&depth, sizeof(uint16_t));  // little-endian
	myNBytes += sizeof(uint16_t);
}

void
swf::
placeObject(uint16_t id, int depth) {
// Replace the object at "depth" with "id"
// PlaceObject2: p 37
	int length = 0;

	// put the body in a separate buffer so we have
	// the length to create a RecordHeader
	char obj[1000];

	// first 8 bits are flags:
	char c = 2;  // PlaceFlagHasCharacter only
	obj[length++] = c;

	// Depth
	memcpy (&obj[length], (char*)&depth, 2);  // little-endian
	length += 2;

	// Id of character to place
	uint16_t idc = id;
	memcpy (&obj[length], (char*)&idc, 2);  // little-endian
	length += 2;

	makeRecordHeader(26, length);
	memcpy(&myBuf[myNBytes], obj, length);
	myNBytes += length;
}

void
swf::
showFrame() {
	makeRecordHeader(1, 0);
}

int
swf::
getCharId() {
	return myId;
}

int
swf::
getNumPixels() {
	return myFramePixels;
}
/////////////////////////// class RecordHeader /////////////////////////////////

RecordHeader::
RecordHeader(int code) :
myCode(code), myLength(-1), myNBytes(-1)
{	}

RecordHeader::
RecordHeader(int code, int length) :
myCode(code), myLength(length)
{	
	if(length > 62 || myCode == 2)
		myNBytes = 6;
	else
		myNBytes = 2;
}

void
RecordHeader::
setLength(int length) {
	myLength = length;

	if (myLength > 62 || myCode == 2)
		myNBytes = 6;
	else
		myNBytes = 2;
}

int
RecordHeader::
write(char * dest) {

	if(myLength < 0 || myNBytes <0){
		cerr << "Error: RecordHeader not fully defined. Cannot write to "
			 << "destination. (Length was never intialized)" << endl;
	}
		
	uint16_t t = myCode << 6;

	// XXX use long tag for code=2 (DefineShape) to match "SWF Uncovered"
	if (myLength > 62 || myCode == 2) {
		t |= 0x3f;
		memcpy(dest, (char*)&t, 2);

		dest += 2;
		int32_t l = myLength;
		memcpy(dest, (char*)&l, 4);

		return myNBytes;

	} else {
		t |= myLength;
		memcpy(dest, (char*)&t, 2);

		return myNBytes;
	}

}

int
RecordHeader::
size() const {
	if(myLength < 0 || myNBytes <0){
		cerr << "Error: RecordHeader has no length or has negative length."
			 << endl;
	}
	return myNBytes;
}
////////////////////////////// class RECT //////////////////////////////////////
Rect::
Rect(int xmin, int xmax, int ymin, int ymax) :
myXmin(xmin), myXmax(xmax), myYmin(ymin), myYmax(ymax)
{
	char nbits = (char)std::max(nsb(xmax), nsb(ymax));
	myNBytes = bytealign(5 + int(4*nbits));
}

// uses copy bits which uses a char *
int
Rect::
write( char *dest) {
	char nbits = (char)std::max(nsb(myXmax), nsb(myYmax));
	int astart = 0;
	int nb = 5;
	int bstart = sizeof(char)*8 - nb;

	copyBits(nb, (char*)&nbits, bstart, dest, astart);
	astart += nb;

	bi((char*)&myXmin, sizeof(myXmin));
	bstart = 8*sizeof(int) - nbits;
	copyBits(nbits, (char*)&myXmin, bstart, dest, astart);
	astart += nbits;

	bi((char*)&myXmax, sizeof(myXmax));
	copyBits(nbits, (char*)&myXmax, bstart, dest, astart);
	astart += nbits;

	bi((char*)&myYmin, sizeof(myYmin));
	copyBits(nbits, (char*)&myYmin, bstart, dest, astart);
	astart += nbits;

	bi((char*)&myYmax, sizeof(myYmax));
	copyBits(nbits, (char*)&myYmax, bstart, dest, astart);
	astart += nbits;

	// return the number of bytes used in buf
	//	myNBytes = bytealign(astart);
	return myNBytes;

}

int
Rect::
size() const {
	return myNBytes;
}
/////////////////////////// class DefineShape //////////////////////////////////

//DefineShape3
DefineShape::
DefineShape( uint16_t shapeId, vector<FlashLine> &lines, 
			vector<FlashPoint> &points) :
  myHeader(32),
  myShapeId(shapeId),
  mySWS(10000,0),
  myEndShape(false),
  myNBytes(0),
  swsNBits(0)
{
	int xmin, xmax, ymin, ymax;
	flashLimits(lines, points, xmin, xmax, ymin, ymax);

	// need to add 20 twips to account for drawing points
	xmax += 20;
	ymax += 20;
	Rect temp(xmin, xmax, ymin, ymax);
	myLimits = temp;

	NumFillBits = 0;
	NumLineBits = 1;

	makeShapeWithStyle(lines, points);

}
void
DefineShape::
makeSCRNewStyles( const int color[4]) {
	// check for resizing
	int currbytes = bytealign(swsNBits);
	if( (mySWS.size() - currbytes) < 100){
		mySWS.resize(mySWS.size()*2, 0);
	}

	char c = 0;

	c = 16;	// StateNewStylesFlag

	copyBits (6, (char*)&c, 2, &mySWS[0], swsNBits);
	swsNBits += 6;
	swsNBits = 8*bytealign(swsNBits);

	makeFillStyleArray();
	makeLineStyleArray(color);
	copyBits(4, (char*) &NumFillBits, 4, &mySWS[0], swsNBits);
	swsNBits += 4;
	copyBits(4, (char*) &NumLineBits, 4, &mySWS[0], swsNBits);
	swsNBits += 4;

} 
//check for resizing in all statements
void
DefineShape::
makeStyleChangeRecord(int x, int y){
	// check for resizing
	int currbytes = bytealign(swsNBits);
	if( (mySWS.size() - currbytes) < 100){
		mySWS.resize(mySWS.size()*2, 0);
	}

	char c = 0;

	c = 9;	// StateLineStyle & StateMoveTo flags on

	copyBits (6, (char*)&c, 2, &mySWS[0], swsNBits);
	swsNBits += 6;

	char moveBits = 14;
	copyBits(5, &moveBits, 3, &mySWS[0], swsNBits);
	swsNBits += 5;

	int16_t moveDeltaX = x;
	bi((char*)&moveDeltaX, sizeof(moveDeltaX));
	copyBits (moveBits, (char*)&moveDeltaX, 16-moveBits, &mySWS[0], swsNBits);
	swsNBits += moveBits;

	int16_t moveDeltaY = y;
	bi((char*)&moveDeltaY, sizeof(moveDeltaY));
	copyBits (moveBits, (char*)&moveDeltaY, 16-moveBits, &mySWS[0], swsNBits);
	swsNBits += moveBits;

	char lsc = 1;
	copyBits(NumLineBits, &lsc, 8-NumLineBits, &mySWS[0], swsNBits);
	swsNBits += NumLineBits;

}

void
DefineShape::
makeEndShapeRecord(){
	int currbytes = bytealign(swsNBits);
	if( (mySWS.size() - currbytes) < 100) {
		mySWS.resize(mySWS.size()*2, 0);
	}

	swsNBits += 6;

	myEndShape = true;
	myHeader.setLength( bytealign(swsNBits) + myLimits.size() + sizeof(myShapeId));
	myNBytes = myHeader.size() + bytealign(swsNBits) + myLimits.size()
				+ sizeof(myShapeId);
}

void
DefineShape::
makeStraightEdgeRecord(int deltaX, int deltaY){
	int currbytes = bytealign(swsNBits);
	if( (mySWS.size() - currbytes) < 100) {
		mySWS.resize(mySWS.size()*2, 0);
	}

	char c = 3;		//TypeFlag=1, StraightFlag = 1
	copyBits(2, &c, 6, &mySWS[0], swsNBits);
	swsNBits += 2;

	int numbits = std::max(nsb(deltaX), nsb(deltaY));		// bits/value
	c = numbits-2;   // number of bits/value = 11 + 2 = 13
	copyBits(4, &c, 4, &mySWS[0], swsNBits);
	swsNBits += 4;

	c = 1;	// general line
	copyBits(1, &c, 7, &mySWS[0], swsNBits);
	swsNBits += 1;

	int16_t dx = deltaX;
	bi((char*)&dx, sizeof(dx));
	copyBits (numbits, (char*)&dx, 16-numbits, &mySWS[0], swsNBits);
	swsNBits += numbits;

	int16_t dy = deltaY;
	bi((char*)&dy, sizeof(dy));
	copyBits (numbits, (char*)&dy, 16-numbits, &mySWS[0], swsNBits);
	swsNBits += numbits;
}

void
DefineShape::
makeVerticalLine( int dy){
	int currbytes = bytealign(swsNBits);
	if( (mySWS.size() - currbytes) < 100) {
		mySWS.resize(mySWS.size()*2, 0);
	}

	char c = 3;		//TypeFlag=1, StraightFlag = 1
	copyBits(2, &c, 6, &mySWS[0], swsNBits);
	swsNBits += 2;

	int numbits = nsb(dy);		// bits/value
	c = numbits-2;   // number of bits/value = 11 + 2 = 13
	copyBits(4, &c, 4, &mySWS[0], swsNBits);
	swsNBits += 4;

	c = 0;	// General line flag ( 0 signals vert/horiz line)
	copyBits(1, &c, 7, &mySWS[0], swsNBits);
	swsNBits += 1;

	c = 1;	// vertical line flag
	copyBits(1, &c, 7, &mySWS[0], swsNBits);
	swsNBits += 1;

	int16_t y = dy;
	bi((char*)&y, sizeof(y));
	copyBits (numbits, (char*)&y, 16-numbits, &mySWS[0], swsNBits);
	swsNBits += numbits;

}

void
DefineShape::
makeHorizLine( int dx){
	int currbytes = bytealign(swsNBits);
	if( (mySWS.size() - currbytes) < 100) {
		mySWS.resize(mySWS.size()*2, 0);
	}
	char c = 3;		//TypeFlag=1, StraightFlag = 1
	copyBits(2, &c, 6, &mySWS[0], swsNBits);
	swsNBits += 2;

	int numbits = nsb(dx);		// bits/value
	c = numbits-2;   // number of bits/value = 11 + 2 = 13
	copyBits(4, &c, 4, &mySWS[0], swsNBits);
	swsNBits += 4;

	c = 0;	// vert/horiz line
	copyBits(1, &c, 7, &mySWS[0], swsNBits);
	swsNBits += 1;

	c = 0;	// horiz line
	copyBits(1, &c, 7, &mySWS[0], swsNBits);
	swsNBits += 1;

	int16_t x = dx;
	bi((char*)&x, sizeof(x));
	copyBits (numbits, (char*)&x, 16-numbits, &mySWS[0], swsNBits);
	swsNBits += numbits;

}

void
DefineShape::
makeRectangleShapeRecord(int width, int height) {

	makeVerticalLine( height);
	makeHorizLine(width);
	makeVerticalLine( -height);
	makeHorizLine( -width);

}

void
DefineShape::
makeShapeWithStyle(const vector<FlashLine> & lines, const vector<FlashPoint> &points){

	int color[4] = {0,0,0,255};
	makeFillStyleArray();
	makeLineStyleArray(color);
	copyBits(4, (char*) &NumFillBits, 4, &mySWS[0], swsNBits);
	swsNBits += 4;
	copyBits(4, (char*) &NumLineBits, 4, &mySWS[0], swsNBits);
	swsNBits += 4;

	makeShapeRecords(lines, points);
}

void
DefineShape::
makeShapeRecords(const vector<FlashLine> & lines, const vector<FlashPoint> &points){

	for(size_t i = 0; i < lines.size(); i++){
		makeSCRNewStyles(lines[i].RGBA0);
		makeStyleChangeRecord(lines[i].x[0], lines[i].y[0]);

		//make delta transforms
		int delX = lines[i].x[1] - lines[i].x[0];
		int delY = lines[i].y[1] - lines[i].y[0];
		makeStraightEdgeRecord(delX, delY);
	}

	for (size_t i = 0 ; i < points.size(); i++){
		makeSCRNewStyles(points[i].RGBA);
		makeStyleChangeRecord(points[i].xC, points[i].yC);
		makeRectangleShapeRecord(20, 20);
	}

	makeEndShapeRecord();
}

void 
DefineShape::
makeFillStyleArray() {
	char FillStyleCount = 0;
	copyBits(8, &FillStyleCount, 0, &mySWS[0], swsNBits);
	swsNBits += 8;
}

void
DefineShape::
makeLineStyleArray(const int color[4]) {
	char LineStyleCount = 1;
	copyBits(8, &LineStyleCount, 0, &mySWS[0], swsNBits);
	swsNBits += 8;
	makeLineStyle(20, color);
	
}

void
DefineShape::
makeLineStyle(uint16_t width, const int RGBA[4]) {

	copyBits(16, (char*)&width, 0, &mySWS[0], swsNBits);
	swsNBits += 16;

	char red = RGBA[0];
	copyBits(8, &red, 0, &mySWS[0], swsNBits);
	swsNBits += 8;

	char green = RGBA[1];
	copyBits(8, &green, 0, &mySWS[0], swsNBits);
	swsNBits += 8;

	char blue = RGBA[2];
	copyBits(8, &blue, 0, &mySWS[0], swsNBits);
	swsNBits += 8;

	char alpha = RGBA[3];
	copyBits(8, &alpha, 0, &mySWS[0], swsNBits);
	swsNBits += 8;

}

int
DefineShape::
write(char* dest) {
	if (!myEndShape) {
		cerr << "attempting to write incomplete shape" <<endl;
	}
		
	dest += myHeader.write(dest);
	// write shape id
	memcpy(dest, &myShapeId, sizeof(myShapeId));
	dest += sizeof(myShapeId);

	dest += myLimits.write(dest);

	//write mySWS for bytealign(swsNBits)
	memcpy(dest, &mySWS[0], bytealign(swsNBits));
	return myNBytes;
}
///////////////////// Additional Functions /////////////////////////////////////
void
parseFeedBuffer(const vector<GLfloat> & buffer, vector<FlashLine>& lines,
			   vector<FlashPoint> &points) {

	points.clear();
	lines.clear();
	for( size_t i=0; i<buffer.size(); i++ ) {
		if( buffer[i] == GL_POINT_TOKEN ) {
			points.push_back( FlashPoint( &buffer[i+1]) );
			i += 7;

		} else if ( buffer[i] == GL_LINE_TOKEN ) {
			lines.push_back( FlashLine( &buffer[i+1]) );
			i += 14;

		} else if ( buffer[i] == GL_LINE_RESET_TOKEN ) {
			lines.push_back( FlashLine( &buffer[i+1]) );
			i += 14;

		} else if (buffer[i] == GL_PASS_THROUGH_TOKEN) {
			break;
		}
	}
}

void
flashLimits(vector<FlashLine>& lines, vector<FlashPoint>& points,
				int& xmin, int& xmax, int& ymin, int& ymax){

	int m = numeric_limits<int>::max();
	xmin = m;
	xmax = -m;
	ymin = m;
	ymax = -m;

	for (size_t i=0; i<lines.size(); i++) {
		xmin = std::min(xmin, lines[i].x[0]);
		xmax = std::max(xmax, lines[i].x[0]);
		xmin = std::min(xmin, lines[i].x[1]);
		xmax = std::max(xmax, lines[i].x[1]);

		ymin = std::min(ymin, lines[i].y[0]);
		ymax = std::max(ymax, lines[i].y[0]);
		ymin = std::min(ymin, lines[i].y[1]);
		ymax = std::max(ymax, lines[i].y[1]);
	}

	for(size_t i=0; i<points.size(); i++) {
		xmin = std::min(xmin, points[i].xC);
		xmax = std::max(xmax, points[i].xC);
		ymin = std::min(ymin, points[i].yC);
		ymax = std::max(ymax, points[i].yC);
	}
}

//move to swf
void
setXYLimits(const vector<GLfloat> &buffer, swf & flashobj) {

	int xmin;
	int xmax;
	int ymin;
	int ymax;
	vector<FlashPoint> points;
	vector<FlashLine> lines;

	parseFeedBuffer(buffer, lines, points);

	flashLimits(lines, points, xmin, xmax, ymin, ymax);

	int scale;
	int framepixels = flashobj.getNumPixels();

	// scale down the image from the openGL feedback buffer to fit
	// within a 300 pixel video
	// NOTE: there are 20 twips in a pixel for a flash video
	if (fabs((double)(xmax-xmin)) > fabs((double)(ymax-ymin))) {
		scale = (framepixels*20)/(xmax-xmin);
	} else {
		scale = (framepixels*20)/(ymax-ymin);
	}

	flashobj.setScale(scale);
	flashobj.setShift(xmin, ymin);

	xmin = 0;
	xmax = scale*(xmax-xmin);
	ymin = 0;
	ymax = scale*(ymax-ymin);
	int framewidth = xmax + xmax/10;
	int frameheight = ymax + ymax/10;

	flashobj.setYMax(ymax);
	flashobj.setFrameDimen(framewidth,frameheight);
}

#ifdef MAIN
#define BUFSIZE 100000
int
main() {
	swf flash;
	int bgcolor[3] = { 100,100,100};
	vector<FlashPoint> somepoint;

	vector<GLfloat> buffer(30);
	buffer[0] = GL_LINE_TOKEN;
	buffer[1] = 100;			//x1
	buffer[2] = 5000; 			//y1
	buffer[3] = 0;				//z1
	buffer[4] = 0;				//R
	buffer[5] = 0;				//G
	buffer[6] = 0;				//B
	buffer[7] = 255;			//A
	buffer[8] = 5000;			//x2
	buffer[9] = 2000;			//y2
	buffer[10] = 0;				//z2
	buffer[11] = 0;				//R
	buffer[12] = 0;				//G
	buffer[13] = 0;				//B
	buffer[14] = 255;			//A
	buffer[15] = GL_LINE_TOKEN;
	buffer[16] = 4000;			//x1
	buffer[17] = 2000;			//y1
	buffer[18] = 0;				//z1
	buffer[19] = 100;			//R
	buffer[20] = 0;				//G
	buffer[21] = 0;				//B
	buffer[22] = 255;			//A
	buffer[23] = 1000;			//x1
	buffer[24] = 400;			//y1
	buffer[25] = 0;				//z1
	buffer[26] = 100;			//R
	buffer[27] = 0;				//G
	buffer[28] = 0;				//B
	buffer[29] = 255;			//A

	flash.setFrameDimen(6000, 6000);
	flash.setBackgroundColor(bgcolor);
	flash.setYMax(6000);
	flash.openBuffer();
	flash.makeShape(buffer);
	flash.placeObject(1,1);
	for(int i =0; i < 40; i++) {
		if (i % 2 == 1)
			flash.removeObject(1);
		if ( i % 2 == 0 && i > 1)
			flash.placeObject(1,1);
		flash.showFrame();
	}
	flash.closeBuffer();
	flash.writeBuffer("test.swf");
}
#endif // MAIN
