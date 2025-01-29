//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
/////////////////////////////////////////////////////////////////////////////
//  Flaps Mode Visualizer flash creation
/////////////////////////////////////////////////////////////////////////////

#ifndef SWF_H
#define SWF_H

#include <vector>
#include <iostream>
#include <stdint.h>
//GLfloat and Tokens are in here
#ifdef __DARWIN__
    #include <OpenGL/gl.h>
#else
    #include <GL/gl.h>
#endif

using namespace std;
/** \class FlashLine
 * \brief Holds information for OpenGL lines to be later used in flash.
 * Ignores z coord.
 */
class FlashLine {
public:
	int x[2], y[2];
	int RGBA0[4];
	int RGBA1[4];

	/** \fn FlashLine()
	 * \brief Default constructor.  Sets every thing to zero.
	 * RGB color is set to black and opaque.
	 */
	FlashLine() {
		 x[0] = y[0] = x[1] = y[1] = 0; 

		 RGBA0[0] = RGBA1[0] = 0;
		 RGBA0[1] = RGBA1[1] = 0;
		 RGBA0[2] = RGBA1[2] = 0;
		 RGBA0[3] = RGBA1[3] = 255;
	}

	/** \fn FlashLine( const GLfloat* buffer)
	 * \brief Create a line from an array/vector.
	 * \param[in] buffer An array of GLfloats that define two vertices.
	 * Should be in the following form: \n
	 * \t x1 y1 z1 R G B A x2 y2 z2 R G B A
	 */
	FlashLine(const GLfloat* buffer ) {
		x[0] = buffer[0];
		y[0] = buffer[1];
		for(size_t i=0; i < 4; i++) {
			RGBA0[i] = buffer[i+3]*255.0;
		}
		cout << endl;
		x[1] = buffer[7];
		y[1] = buffer[8];
		for(size_t i=0; i < 4; i++) {
			RGBA1[i] = buffer[i+10]*255.0;
		}
	}
};

/** \class FlashPoint
 * \brief Creates a point for flash video.  Ignores z coord.
 */
class FlashPoint{
 public:
	int xC,yC;
	int RGBA[4];

	/** \fn FlashPoint()
	 * \brief Default Constructor. Set (x,y) = (0,0) and color to black
	 */
	FlashPoint(){ 
		xC = yC = 0; 
		RGBA[0] = 0;
		RGBA[1] = 0;
		RGBA[2] = 0;
		RGBA[3] = 255;
	}

	/** \fn FlashPoint(const GLfloat* buf)
	 * \brief Create a point from an array of GLfloats generated from 
	 * feedback buffer mode in OpenGL.
	 * \param[in] buf array of GLfloats that defines a vertex. Order should
	 * be as follows: \n
	 * \t x y z R G B A
	 */
	FlashPoint( const GLfloat* buf ) {
		xC = (int) buf[0];
		yC = (int) buf[1];

		for(size_t i=0; i < 4; i++) {
			RGBA[i] = buf[i+3]*255.0;
		}
		
	}
};

using std::vector;
using std::string;

// makes a flash video from an openGL feedbuffer
// scales the openGL coords and performs transformation since openGL and
// flash have different coordinate systems
class swf {
	private:
		vector<char> myBuf;					// Buffer that hold flash code
		uint32_t myNBytes;					// Total buffer size
		uint16_t myId;						// Shape Id counter

		uint16_t myTime;					// Video time
		uint16_t myNframes;					// Number of frames
		uint16_t myFps;						// Frames per second

		int myScale;						// scale factor between flash and openGL
		int myXShift;						// Shifting left
		int myYShift;						// Shifting up
		int myYMax;							// for inverting OpenGl coords (Max over all frames)
		int myFrameWidth;					// width of video in twips
		int myFrameHeight;					// height of video in twips
		int myFramePixels;					// max number of pixels for video
		int myBGColor[3];					// background color

		// private member functions
		void makeRecordHeader(int code, int length);
		void makeBackgroundColor(int color[]);
		void makeFileHeader(int xmin, int xmax, int ymin, int ymax, 
							uint16_t nframes);
		
		void makeRect(int xmin, int xmax, int ymain, int ymax);
		void changeCoords(vector<FlashPoint> & points);
		void changeCoords(vector<FlashLine> & lines);

	public:
		// constructor
		swf();

		// Modifier

		/** \fn void setScale(int scaleVal)
		 * \brief Set the Scale.
		 */
		void setScale(int scaleVal);

		/** \fn void setBackgroundColor(int color[])
		 * \brief Sets the background color to the RGB values in color.
		 * \details Default background color is white. Background color
		 * must be set prior to opening the buffer.  Once the buffer is opened
		 * changing the background color does nothing.
		 */
		void setBackgroundColor(int color[]);

		/** \fn void setNframes(int numFrames)
		 * \brief Change the number of frames in the flash file.
		 * \details Default number of frames is 40.  This function will
		 * also change the frame rate.  It is the user's responsibility
		 * to choose the number of frames 
		 */
		void setNframes(int numFrames);

		/** \fn void setFramePixels(int pixelCount)
		 * \brief Changes the maximum number of pixels that the video can be.
		 * \details Increasing the frame pixels will likely lead to larger video
		 * sizes and different scaling.  Frame pixels are used in determining
		 * the frame height, frame width, and scaling.
		 */
		void setFramePixels(int pixelCount);

		/** \fn void setFrameDimen(int frameWidth, int frameHeight)
		 * \brief Set the width and height of the video.
		 * \details Numbers are specified in twips.  These are variables that must
		 * be defined prior to opening the buffer.
		 */
		void setFrameDimen(int frameWidth, int frameHeight);

		/** \fn void setTime(int seconds)
		 * \brief Set length of video in seconds.
		 */
		void setTime(int seconds);

		/** \fn void setShift(int XShift, int YShift)
		 * \brief Units of shift in X and Y.  Specified in OpenGL coords.
		 * \details Shifting is done in order to get the scaling correct and
		 * decrease white space. 
		 */
		void setShift(int XShift, int YShift);

		/** \fn void setYMax(int YMax)
		 * \brief Set the maximum y value encountered across all frames.
		 * \details YMax is used to invert the buffer from OpenGL.  This needs
		 * to be done because the origin in OpenGL is different from that of
		 * Flash.
		 */
		void setYMax(int YMax);

		/** \fn void openBuffer()
		 * \brief Opens the buffer for writing
		 * \details Certain variables must be set prior to opening the buffer.
		 * It is the user's responsibility to make sure that these values are
		 * defined.  No warnings or errors will be given. \n \n
		 * The variables are as follows: \n
		 * background color, framewidth, frameheight, number of frames, frames 
		 * per second.
		 */
		void openBuffer();

		/** \fn void closeBuffer()
		 * \brief Close the buffer for writing.  Writes the file size back into
		 * the file header.
		 */
		void closeBuffer();

		/** \fn void writeBuffer(string path)
		 * \brief Write the buffer out to the filename path.
		 * \details Uses memcpy.
		 */
		void writeBuffer(string path);

		/** \fn void makeShape( const vector<GLfloat> & buffer)
		 * \brief Create a shape using the buffer from OpenGL feedback mode.
		 * \details Uses a DefineShape3.  It creates a shape writes it to
		 * myBuf.
		 * \param[in] buffer The end of the buffer should be marked with a 
		 * glPassThrough to allow parsefeedbuffer to determine the end.
		 */
		void makeShape( const vector<GLfloat> & buffer);

		/** \fn void removeObject(int depth)
		 * \brief Removes all placed objects at specified depth.
		 * \details This is RemoveObject2 in the SWF Manual
		 */
		void removeObject(int depth);

		/** \fn void placeObject(uint16_t id, int depth)
		 * \brief Place an object with "id" at "depth"
		 * \details Should be PlaceObject2 in swf manual
		 */
		void placeObject(uint16_t id, int depth);

		/** \fn void showFrame()
		 * \brief Same as ShowFrame which is described in the manual.
		 */
		void showFrame();

		//accessors
		/** \fn int getCharId()
		 * \brief Returns the id of the next shape to be created.
		 */
		int getCharId();
		
		/** \fn int getNumPixels()
		 * \brief Return myFramePixels.
		 */
		int getNumPixels();

};

/** \class RecordHeader
 * \brief Creates a RecordHeader as noted in swf manual.
 * \details The class will automatically determine the size and format
 * based on the length.
 */
class RecordHeader{
	private:
		int myCode;
		int myLength;
		int myNBytes;

	public:
		/** \fn RecordHeader( int code)
		 * \brief Constructor: Used when length is not yet determined.
		 * \details Primary use is for the class DefineShape.
		 */
		RecordHeader(int code);

		/** \fn RecordHeader(int code, int length)
		 * \brief Constructor: Creates a RecordHeader when both the code and
		 * length are known.
		 */
		RecordHeader(int code, int length);

		/** \fn void setLength(int length)
		 * \brief Set the length of the RecordHeader if it was not previously set.
		 */
		void setLength(int length);

		/** \fn int write(char *dest)
		 * \brief Write the RecordHeader to an array specified by dest.
		 * \return Number of bytes writen to destination
		 */
		int write(char *dest);

		/** \fn int size() const
		 * \brief Size of the RecordHeader
		 * \return Returns either 2 bytes or 6 as specified by the SWF manual.
		 */
		int size() const;
};

/** \class Rect
 * \brief Creates a RECT that is used to describe the limits in cartesian coordinate
 * systems.
 */
class Rect {
	private:
		int myXmin, myXmax, myYmin, myYmax;
		int myNBytes;

	public:
		/** \fn RECT()
		 * \brief Constructor: Default Constructor.
		 * \details Leave everything unitialized.
		 */
		Rect() {}

		/** \fn RECT(int xmin, int xmax, int ymin, int ymax)
		 * \brief Creates a RECT with the the limits specified.
		 */
		Rect(int xmin, int xmax, int ymin, int ymax);

		/** \fn int write(char *dest)
		 * \brief Writes a RECT to a character array dest.
		 * \pre dest should have enough memory allocated such that RECT
		 * can be written without segmentation fault.
		 * \returns Number of bytes written to destination
		 */
		int write(char *dest);

		/** \fn int size() const
		 * \brief Returns the number of bytes required to write RECT.
		 */ 
		int size() const;

};

/** \class DefineShape
 * \brief Creates a DefineShape3 from FlashLine and FlashPoint vectors.
 * \details Each define shape is used to represent a single phase angle
 * as the video cycles from 0 degrees phase all the way to 360. Points 
 * in the DefineShape are actually represented by four individual line
 * segments.
 */
class DefineShape {
	private:
		RecordHeader myHeader;
		uint16_t myShapeId;		// ID for flash dictionary
		Rect myLimits;			// limits for the shape
		vector<char> mySWS;		// SWS = SHAPEWITHSYTLE
		bool myEndShape;		// Boolean for if EndShapeRecord has been written
		int myNBytes;			// number of total bytes 
		int swsNBits;			// Number of bits written to mySWS
		int NumLineBits;		// Number of Line bits for LINESTYLEARRAY
		int NumFillBits;		// Number of fill bits for FILLSYTLEARRAY

		// private functions

		/** \fn void makeSCRNewStyles( const int color[4])
		 * \brief Creates new fill and line style arrays only.
		 * \details The function creates new fill and line styles
		 * in the shape record array. SCR stands for Style Change Record.
		 * \param[in] color RGBA values between 0 and 255
		 */
		void makeSCRNewStyles( const int color[4]);

		/** \fn void makeStyleChangeRecord( int x, int y)
		 * \brief Move to the specified (x,y) coordinate.
		 * \details The StyleChangeRecord that gets created with move to the
		 * Flash cartesian point (x,y) and select a line style from the line
		 * style array.
		 */
		void makeStyleChangeRecord( int x, int y);

		/** \fn void makeEndShapeRecord()
		 * \brief Ends the shape record array.
		 * \details Makes the assumption that mySWS was intialized to zero.
		 */
		void makeEndShapeRecord();

		/** \fn void makeStraightEdgeRecord(int deltaX, int deltaY)
		 * \brief Draws a line from current point with the deltas
		 * \details Positive deltaX is to the right and positive
		 * deltaY is down. Used to draw the lines in the flash file
		 */
		void makeStraightEdgeRecord(int deltaX, int deltaY);

		/** \fn void makeVerticalLine(int dy)
		 * \brief Make a vertical line
		 * \details positve dy is down since origin is in the top left
		 */
		void makeVerticalLine(int dy);

		/** \fn void makeHorizLine(int dx)
		 * \brief Create a horizontal line shape record.
		 * \details positive dx is to the right since origin is in the
		 * top left.
		 */
		void makeHorizLine(int dx);

		/** \fn void makeRectangleShapeRecord(int width, int height)
		 * \brief Create a rectangle of "width" and "height".
		 * \details Rectangle is created down and to the right of the current
		 * point. This function uses the horizontal and vertical line functions
		 * and is used to draw points.
		 */
		void makeRectangleShapeRecord(int width, int height);

		/** \fn void makeShapeWithStyle(const vector<FlashLine> & lines,
									const vector<FlashPoint> & points)
		 * \brief Creates a Shape with style.
		 * \detail This function creates a ShapeWithStyle that has the
		 * NumFillBits, NumLineBits, FILLSTYLEARRAY, LINESTYLEARRAY, and
		 * ShapeRecords. Called exclusively from the DefineShape constructor.
		 * \sa makeShapeRecords()
		 * \param[in] lines Vector of lines that are in the DefineShape
		 * \param[in] points Vector of points that are in the DefineShape.
		 */
		void makeShapeWithStyle(const vector<FlashLine> & lines,
								const vector<FlashPoint> & points);

		/** \fn void makeShapeRecords(const vector<FlashLine>& lines, 
		 *							  const vector<FlashPoint>& points)
		 * \brief Creates the ShapeRecords for makeShapeWithStyle.
		 * \details Examples of shape records include Style Change Records, 
		 * End Shape Record, Straight Line Records, etc. Called exclusively
		 * from make ShapeWithStyle.
		 * \param[in] lines vector of lines for the shape records.
		 * \param[in] points vector of points for the shape record.
		 */
		void makeShapeRecords(const vector<FlashLine>& lines,
							  const vector<FlashPoint> & points);

		/** \fn void makeFillStyleArray()
		 * \brief Creates a FillStyleArray
		 * \details For the purposes of this flash file there are no solidly
		 * filled shapes so the fill style array is empty.
		 */
		void makeFillStyleArray();

		/** \fn void makeLineStyleArray(const int color[4])
		 * \brief Creates a LineStyleArray.
		 * \details Make sure to convert RGBA numbers from OpenGL. This
		 * is currently being done in FlashPoint and FlashLine constructors.\n
		 * This function is called by both the makeShapeWithStyle and 
		 * makeSCRNewStyles
		 *\param[in] color RGBA numbers between 0 and 255
		 */
		void makeLineStyleArray(const int color[4]);

		/** \void makeLineStyle(uint16_t width, const int RGBA[4])
		 * \brief Creates a single line style for a LineStyleArray.
		 * \param[in] width Width of the lin in twips
		 * \param[in] RGBA color and alpha value for the line
		 */
		void makeLineStyle(uint16_t width, const int RGBA[4]);

	public:

		/** \fn DefineShape(uint16_t shapeId, vector<FlashLine> &lines,
		 *					vector<FlashPoint> &points)
		 * \brief
		 */
		DefineShape(uint16_t shapeId, vector<FlashLine> &lines, 
					vector<FlashPoint> & points);

		/** \fn int write(char* dest)
		 * \brief Write a define shape to the specified character array.
		 * \pre Array should be large enough to allow the DefineShape to be
		 * written without causing segmentation faults.
		 */
		int write(char *dest);

		/** \fn int size() const
		 * \brief Returns the total number of bytes needed to write
		 * the DefineShape.
		 * \details This number represents the combination of the RecordHeader,
		 * shape id, rect, and shape with style.
		 */
		int size() const { return myNBytes; }
};

/** \fn void parseFeedBuffer(const vector<GLfloat> & buffer, vector<FlashLine> & lines,
			   vector<FlashPoint> &points)
  * \brief Takes a vector from GL_Feedback_buffer mode and creates FlashLines and 
  * FlashPoints from it.
  */
void 
parseFeedBuffer(const vector<GLfloat> & buffer, vector<FlashLine> & lines,
			   vector<FlashPoint> &points);

/** \void flashLimits(vector<FlashLine>& lines, vector<FlashPoint>& points,
 *					int& xmin, int& xmax, int& ymin, int& ymax);
 * \brief Given vectors of lines and points find the max and min cartesian
 * coordinates in the set.
 */
void
flashLimits(vector<FlashLine>& lines, vector<FlashPoint>& points,
			int& xmin, int& xmax, int& ymin, int& ymax);

/** fn void setXYLimts(const vector<GLfloat>& buffer, swf &flashobj)
 * \brief Sets the scale and applicable limits for a flash object given
 * \detail The scale, frame dimensions, shifting, and ymax for a flash
 * object are all determined and set in this function.  The buffer should 
 * contain phase angles where all the mode shape is at its max. Ie 0, 90, 180
 * 270.
 * \note should consider moving inside of the class swf
 * \param[in] buffer
 * \param[in] flashobj 
 */
void
setXYLimits(const vector<GLfloat> &buffer, swf & flashobj);

// Use DEFINESHAPE3
// How do I add colors?  STYLECHANGERECORD Sufficient?

#endif // SWF_H

