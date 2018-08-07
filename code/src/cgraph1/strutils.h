#ifndef STRUTILS_H_INCLUDED
#define STRUTILS_H_INCLUDED

#include <string.h>
#include <ctype.h>

/**
 * removes spaces at the beginning
 * of the string
 **/
char *str_remove_sps_sol( char *dest, const char *str );

/**
 * receives a fileName with path
 * and returns only the fileName
 **/
char *getFileName(char *destiny, const char *fileWithPath);

char* getParamName(char* target, const char* str);

char* getParamValue(char* target, const char* str);

/* converts all chacarters to uppercase */
void str_all_to_upper( char *str );

/* removes spaces and other meaningless characthers from end of line */
char* str_remove_sps_eol( char *str );

/* char is invisible:  if char = ' ' or '\n' or '\t' or '\r' */
char char_is_invisible( const char c );

/* to align at right */
void str_fill_spaces_left( char *dest, const char *str, int n );

/* to align at left */
void str_fill_spaces_right( char *dest, const char *str, int n );

/* to align at center */
void str_fill_spaces_both( char *dest, const char *str, int n );

/* remove repeated spacing characters (e.g. spaces of tabs ) */
char *str_remove_dbl_spaces( char *str );

char *str_clear_spaces(char *str);

/**
 * splits a string usign a delimiter char
 * returns how many columns were found
 * multDel = 1 indicates that multiple delimiters
 * grouped together represent only one
 **/
int splitString( char **columns, const char *str, const char delimiter,
                 const int maxColumns, const int columnSize, const char multDel );

/**
 * number of numeric characters in
 * line
 **/
int digitsInLine( const char *str, int strSize );


#endif
