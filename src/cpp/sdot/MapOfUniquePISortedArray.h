#pragma once

#include <vfs/support/TypeConfig.h>
#include <map>

/**

*/
template<class T,int s,class Val,bool s_is_neg = ( s < 0 )>
struct MapOfUniquePISortedArray;

///
template<class T,int s,class Val>
struct MapOfUniquePISortedArray<T,s,Val,true> {
};

///
template<class T,int s,class Val>
struct MapOfUniquePISortedArray<T,s,Val,false> {
    using Vals           = std::map<std::array<T,s>,Val>;

    void  init           ( Vfs::PI max_PI_value, const Val &default_value = {}) { this->default_value = default_value; }

    Val&  operator[]     ( std::array<T,s> a ) { auto iter = values.find( a ); if ( iter == values.end() ) iter = values.insert( iter, { a, default_value } ); return iter->second; }

    Val   default_value; ///<
    Vals  values;        ///<
};

///
template<class T,class Val>
struct MapOfUniquePISortedArray<T,0,Val> {
    void init             ( Vfs::PI /*max_PI_value*/, const Val &default_value = {} ) { value = default_value; }

    Val& operator[]       ( std::array<T,0> ) { return value; }

    Val  value;           ///<
};

///
template<class T,class Val>
struct MapOfUniquePISortedArray<T,1,Val> {
    using Vals            = std::vector<Val>;

    void  init            ( Vfs::PI max_PI_value, const Val &default_value = {} ) { values.resize( max_PI_value ); for( Val &val : values ) val = default_value; }

    Val&  operator[]      ( std::array<T,1> a ) { return values[ a[ 0 ] ]; }

    Vals  values;          ///<
};

///
template<class T,class Val>
struct MapOfUniquePISortedArray<T,2,Val> {
    using Vals            = std::vector<Val>;

    void  init            ( Vfs::PI max_PI_value, const Val &default_value = {} ) { values.resize( ( max_PI_value - 1 ) * max_PI_value / 2 ); for( Val &val : values ) val = default_value; }

    Val&  operator[]      ( std::array<T,2> a ) { return values[ ( a[ 1 ] - 1 ) * a[ 1 ] / 2 + a[ 0 ] ]; }

    Vals  values;          ///<
};
