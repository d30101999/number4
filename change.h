#ifndef CHANGE_H
#define CHANGE_H

#include <iostream>

void write_header(FILE *file_out, char char_header, int width, int height, unsigned int max_value);
void write_data(FILE *file_out, int k_bytes, unsigned char *pix_data);
int read_header(char &char_header, int &width, int &height, unsigned int &max_value, FILE *file_in);
void read_data(FILE *file_in, int k_bytes, unsigned char *pix_data);
void free_data(FILE *file_in, FILE *file_out, unsigned char *pix_data);
void free_data(FILE *file, unsigned char *pix_data);
void free_data(FILE *file);
void free_data(unsigned char *pix_data);
namespace color {
    void write_header(FILE *file_out, int width, int height, unsigned int max_value);
    void write_to_file(FILE *file_out, int width, int height, unsigned int max_value,
                       unsigned char *pix_data);
}

namespace gray {
    void write_header(FILE *file_out, int width, int height, unsigned int max_value);
    void write_to_file(FILE *file_out, int width, int height, unsigned int max_value,
                       unsigned char *pix_data);
}

void write_to_file(FILE *file_out, char char_header, int width, int height, unsigned int max_value,
                   unsigned char *pix_data);

unsigned char limit_brightness(double brightness);

double to_sRGB(double _brightness);
double to_sRGB(int brightness);
double from_sRGB(double _brightness);
double from_sRGB(int brightness);

double change_pix_gamma_to_print(double _brightness, double gamma);
double change_pix_gamma_to_print(unsigned char pix_data, double gamma);

double change_pix_gamma_from_file(double _brightness, double gamma);
double change_pix_gamma_from_file(unsigned char pix_data, double gamma);

void draw_pix(unsigned char *pix_data, int width, int x, int y, int brightness, double gamma);
void draw_pix(unsigned char *pix_data, int width, int x, int y, double _brightness, double gamma);

void decode_gamma_from_file(unsigned char *pix_data, int k, double gamma);

#endif // CHANGE_H
