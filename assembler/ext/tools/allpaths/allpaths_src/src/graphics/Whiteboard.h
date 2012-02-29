#ifndef WHITEBOARD_H
#define WHITEBOARD_H

#include "graphics/Color.h"
#include "graphics/Eps.h"

#include "String.h"
#include "Vec.h"

namespace ns_whiteboard {

  typedef pair<float,float> xy_coords;

  enum vert_align { align_top, align_middle, align_bottom };
  enum horiz_align { align_left, align_center, align_right };

class display_type
{
 public:
  virtual void DrawPoint( const xy_coords &coord ) = 0;
  virtual void DrawPoint( double x, double y ) = 0;
  virtual void DrawLine( const xy_coords &startCoord,
			 const xy_coords &stopCoord ) = 0;
  virtual void DrawRect( const xy_coords &topLeftCoord,
			 const xy_coords &bottomRightCoord ) = 0;
  virtual void DrawText( const xy_coords& coords,
                         const String& chars,
                         float rot_angle,
                         float point_size,
                         const vert_align vAlign,
                         const horiz_align hAlign ) = 0;
  virtual void DrawArc( const xy_coords &coords,
			const float radius,
			const float startAngle,
			const float stopAngle ) = 0;

  virtual const color& GetColor() const = 0;
  virtual void SetColor( const color& C ) = 0;

  virtual float GetLineWidth() const = 0;
  virtual void SetLineWidth( const float width ) = 0;

  virtual const String& GetFontName() const = 0;
  virtual float         GetFontSize() const = 0;
  
  virtual void SetFont( const String& font, const float size ) = 0;

  virtual void Finish() = 0;
public:
  virtual ~display_type(){}
};


class graphic 
{
 public:
  graphic() {}

  virtual void Draw( display_type *d ) = 0;
public:
  virtual ~graphic(){}
};


class whiteboard
{
 public:
  whiteboard( ) {}
 
  void Add( graphic *g )
    {
      m_graphics.push_back(g);
    }

  void DisplayOn( display_type *d )
    {
      for ( int i=0; i<(int) m_graphics.size(); ++i )
	m_graphics[i]->Draw( d );
      d->Finish();
    }

  ///Clear all data, deleting all added pointers.
  void DeletePointers() {
    for (int i=0; i != m_graphics.isize(); ++i) {
      delete m_graphics[i];
    }
    m_graphics.clear();
  }

  /// Number of graphics pointers contained.
  int size() const { return m_graphics.size(); }

  ///Clear all data but do not delete pointers: creator is responsible 
  /// for them.
  void clear() { m_graphics.clear(); }

 private:

  vec<graphic *> m_graphics;
};



class point : public graphic
{
 public:
  point( const xy_coords& coords,
	 const float size = 1.0,
	 const color& c = black ) 
    : m_coords(coords),
      m_size(size),
      m_color(c) {}

  void Draw( display_type *d )
    {
      d->SetLineWidth(m_size);
      d->SetColor(m_color);
      d->DrawPoint(m_coords);
    }

 protected:
  xy_coords m_coords;
  float m_size;
  color m_color;
public:
  virtual ~point(){}
};

///Point class with double coordinates to avoid rounding errors.
class dpoint : public graphic
{
 public:
  dpoint( double x, double y,
	 const float size = 1.0,
	 const color& c = black ) 
    : x_(x), y_(y),
      m_size(size),
      m_color(c) {}

  void Draw( display_type *d )
    {
      d->SetLineWidth(m_size);
      d->SetColor(m_color);
      d->DrawPoint(x_,y_);
    }

 protected:
  double x_;
  double y_;
  float m_size;
  color m_color;
public:
  virtual ~dpoint(){}
};


class line : public graphic
{
 public:
  line( const xy_coords& startCoords, 
	const xy_coords& stopCoords,
	const float width = 1.0,
	const color& c = black ) 
    : m_startCoords(startCoords),
      m_stopCoords(stopCoords),
      m_width(width),
      m_color(c)
    {}

  void Draw( display_type *d )
    {
      d->SetLineWidth(m_width);
      d->SetColor(m_color);
      d->DrawLine(m_startCoords,m_stopCoords);
    }

 protected:
  xy_coords m_startCoords, m_stopCoords;
  float m_width;
  color m_color;
public:
  virtual ~line(){}
};



class rect : public graphic
{
 public:
  rect( const xy_coords& topLeftCoords, 
	const xy_coords& bottomRightCoords,
	const color& c = black ) 
    : m_topLeftCoords(topLeftCoords),
      m_bottomRightCoords(bottomRightCoords),
      m_color(c)
    {}

  void Draw( display_type *d )
    {
      d->SetColor(m_color);
      d->DrawRect(m_topLeftCoords,m_bottomRightCoords);
    }

 protected:
  xy_coords m_topLeftCoords, m_bottomRightCoords;
  color m_color;
public:
  virtual ~rect(){}
};



// coords : x-y coordinates of the arc's center of curvature
// radius : radius of curvature
// startAngle, stopAngle : measured counterclockwise from positive x-axis
class arc : public graphic
{
 public:
  arc( const xy_coords & coords,
       const float radius,
       const float startAngle,
       const float stopAngle,
       const float width = 1.0,
       const color &c = black )
    : m_coords(coords),
      m_radius(radius),
      m_startAngle(startAngle),
      m_stopAngle(stopAngle),
      m_width(width),
      m_color(c)
    {}

  void Draw( display_type *d )
    {
      d->SetLineWidth(m_width);
      d->SetColor(m_color);
      d->DrawArc(m_coords,m_radius,m_startAngle,m_stopAngle);
    }

 protected:
  xy_coords m_coords;
  float m_radius, m_startAngle, m_stopAngle, m_width;
  color m_color;
public:
  virtual ~arc(){}
};


class text : public graphic
{
 public:
  text( const xy_coords& coords,
	const String& chars,
	const color& c = black,
	const float fontSize = 12.0,
	const String& fontName = "Times-Roman",
	const float rotAngle = 0.0 )
    : m_coords(coords),
      m_chars(chars),
      m_color(c),
      m_fontSize(fontSize),
      m_fontName(fontName),
      m_rotAngle(rotAngle),
      m_vertAlign(align_bottom),
      m_horizAlign(align_center)
  {}

  void SetVertAlign( const vert_align a ) { m_vertAlign = a; }
  void SetHorizAlign( const horiz_align a ) { m_horizAlign = a; }

  virtual ~text() {}

  void Draw( display_type *d )
    {
      d->SetFont(m_fontName,m_fontSize);
      d->SetColor(m_color);
      d->DrawText(m_coords,m_chars,m_rotAngle,m_fontSize,m_vertAlign,m_horizAlign);
    }

 protected:
  xy_coords m_coords;
  String m_chars;
  color m_color;
  float m_fontSize;
  String m_fontName;
  float m_rotAngle;
  vert_align m_vertAlign;
  horiz_align m_horizAlign;
};




class decorator : public graphic
{
 public:
  decorator(graphic *g) : mp_graphic(g) {}
  
  virtual void Draw( display_type *d )
    { 
      mp_graphic->Draw(d);
    }

 protected:
  graphic *mp_graphic;

public:
  virtual ~decorator(){}
};



class ps_display : public display_type
{
 public:
  ps_display( ostream &ostrm,
	      double horizSize, double vertSize,
	      float border = 0.0 )
    : m_outstrm(ostrm),
      m_color(white),
      m_lineWidth(0.0)
    {
      PrintEpsHeader( m_outstrm, horizSize, vertSize, border ); 
    }

  virtual ~ps_display() {}
  
  void DrawPoint( const xy_coords& coords)
    {
      m_outstrm <<coords.first-m_lineWidth/2<<" "<<coords.second<<" moveto" << "\n";
      m_outstrm <<m_lineWidth<<" 0 rlineto"<< "\n";
      m_outstrm <<"stroke" << "\n";
    }

  void DrawPoint( double x, double y)
    {
      m_outstrm << x - m_lineWidth/2<<" "<< y <<" moveto" << "\n";
      m_outstrm << m_lineWidth << " 0 rlineto"<< "\n";
      m_outstrm << "stroke" << "\n";
    }

  void DrawLine( const xy_coords& startCoords,
		 const xy_coords& stopCoords )
    {
      m_outstrm <<startCoords.first<<" "<<startCoords.second<<" moveto" << "\n";
      m_outstrm <<stopCoords.first<<" "<<stopCoords.second<<" lineto" << "\n";
      m_outstrm <<"stroke" << "\n";
    }

  void DrawText( const xy_coords& coords,
		 const String& chars,
		 float rot_angle,
                 float point_size,
                 const vert_align vAlign,
                 const horiz_align hAlign )
    {
      m_outstrm <<"("<<chars<<") newpath" << "\n";
      m_outstrm <<coords.first<<" "<<coords.second<<" moveto" << "\n";

      if( rot_angle != 0 )
	m_outstrm << "gsave " << rot_angle << " rotate ";

      switch( hAlign ) {
        case align_left:
          // do nothing
          break;
        case align_center:
          m_outstrm <<"("<<chars<<") stringwidth pop 2 div neg 0 rmoveto ";
          break;
        case align_right:
          m_outstrm <<"("<<chars<<") stringwidth pop neg 0 rmoveto ";
          break;
      }

      switch ( vAlign ) {
        case align_top:
          // move down 80% of the font size
          m_outstrm << "0 " << point_size << " 0.8 mul neg rmoveto ";
          break;
        case align_middle:
          // move down 40% of the font size
          m_outstrm << "0 " << point_size << " 0.4 mul neg rmoveto ";
          break;
        case align_bottom:
          // do nothing
          break;
      }

      m_outstrm << "show ";

      if( rot_angle != 0 )
	m_outstrm << " grestore";

      m_outstrm << "\n";
    }

  void DrawArc( const xy_coords &coords,
		const float radius,
		const float startAngle,
		const float stopAngle )
    {
      m_outstrm <<"newpath" << "\n";
      m_outstrm << coords.first <<" "<<coords.second<<" "<<radius<<" "<<startAngle<<" "
		<< stopAngle<<" arc"<< "\n";
      m_outstrm <<"stroke" << "\n";
    }

  void DrawRect( const xy_coords &topLeftCoords,
                 const xy_coords &bottomRightCoords )
    {
      m_outstrm <<"newpath" << "\n";
      m_outstrm << topLeftCoords.first <<" "<< topLeftCoords.second << " moveto" << "\n";
      m_outstrm << bottomRightCoords.first <<" "<< topLeftCoords.second << " lineto" << "\n";
      m_outstrm << bottomRightCoords.first <<" "<< bottomRightCoords.second << " lineto" << "\n";
      m_outstrm << topLeftCoords.first <<" "<< bottomRightCoords.second << " lineto" << "\n";
      m_outstrm << topLeftCoords.first <<" "<< topLeftCoords.second << " lineto" << "\n";
      m_outstrm <<"fill" << "\n";
    }
		

  const color& GetColor() const { return m_color; }

  void SetColor( const color& C )
    {
      if ( C != m_color )
	{
	  m_outstrm <<C.R()<<" "<<C.G()<<" "<<C.B()<<" setrgbcolor" << "\n";
	  m_color = C;
	}
    }

  float GetLineWidth() const { return m_lineWidth; }

  void SetLineWidth( const float width )
    {
      if ( width != m_lineWidth )
	{
	  m_outstrm << width << " setlinewidth" << "\n";
	  m_lineWidth = width;
	}
    }

  const String& GetFontName() const { return m_fontName; }
  float         GetFontSize() const { return m_fontSize; }
  
  void SetFont( const String& font, const float size )
    {
      if ( font != m_fontName || size != m_fontSize )
	{
	  m_outstrm << "/"<<font<<" findfont "<<size<<" scalefont setfont" << "\n";
	  m_fontName = font;
	  m_fontSize = size;
	}
    }

  void Finish() { m_outstrm << flush; }

 protected:
  ostream& m_outstrm;
  color m_color;
  float m_lineWidth;
  String m_fontName;
  float m_fontSize;
};

}

#endif
