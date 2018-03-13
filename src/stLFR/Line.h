#ifndef __STLFR_LINE_H__
#define __STLFR_LINE_H__

#include <vector>
#include <functional>
namespace BGIQD {
namespace stLFR {
typedef std::vector<unsigned int>  base_line;

//float ScoreBaseLine( const base_line & line, 

class ILine
{
    public:
        typedef std::function<float ( const base_line & )> LineFactor;
        virtual base_line LineSelect( LineFactor func ) = 0;

    public :
        virtual ~ILine() {}
};


class Line : public ILine
{
    public:
        virtual base_line LineSelect(LineFactor)
        {
            return m_line;
        }

    public :
        virtual ~Line(){}
    private :
        base_line  m_line;
};

class MultiLine : public ILine
{
    public:
        virtual base_line LineSelect(LineFactor func)
        {
            if(! sorted )
            {
                SelectBySort(func);
                sorted = true ;
            }
            return m_line;
        }

    public :
        virtual ~MultiLine()
        {
            for( int i = 0 ; i<(int) m_lines.size(); i++)
            {
                if( m_lines[i] != nullptr )
                {
                    delete m_lines[i];
                    m_lines[i] = nullptr;
                }
            }
        }
    private :
        void SelectBySort( LineFactor & factor )
        {
            int index = 0;
            float highest = 0.0f ;
            for( int i = 0 ; i< (int) m_lines.size() ; i++)
            {
                float ret = factor(m_lines[i]->LineSelect(factor));

                if( ret >= highest )
                {
                    highest = ret ;
                    index = i;
                }
            }
            m_line = (m_lines[index])->LineSelect(factor);
        }
    private:
        std::vector<Line*> m_lines;
        base_line  m_line;
        bool sorted;
}; //class MultiLine;


class LineFactory
{
    public:

};


}// namespace stLFR
}// namespace BGIQD
#endif //__STLFR_LINE_H__
