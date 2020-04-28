#include "change.h"
#include <cmath>

void write_header(FILE *file_out, char char_header, int width, int height, unsigned int max_value) {
    fseek(file_out, 0, SEEK_SET);
    fprintf(file_out, "P%c\n%i %i\n%i\n", char_header, width, height, max_value);
}

void write_data(FILE *file_out, int k_bytes, unsigned char *pix_data) {
    fwrite(pix_data, k_bytes, 1, file_out);
}

int read_header(char &char_header, int &width, int &height, unsigned int &max_value, FILE *file_in) {
    int scanned = fscanf(file_in, "P%c\n%i %i\n%i\n", &char_header, &width, &height, &max_value);
    if (scanned != 4)
        return 1;
    if (width <= 0 || height <= 0 || max_value <= 0 || max_value > 255 || (char_header != '5' && char_header != '6'))
        return 2;
    return 0;
}

void read_data(FILE *file_in, int k_bytes, unsigned char *pix_data) {
    if (!pix_data)
        return;
    fread(pix_data, 1, k_bytes, file_in);
}

void free_data(FILE *file_in, FILE *file_out, unsigned char *pix_data) {
    if (file_in) fclose(file_in);
    if (file_out) fclose(file_out);
    if (pix_data) free(pix_data);
}

void free_data(FILE *file, unsigned char *pix_data) {
    free_data(nullptr, file, pix_data);
}

void free_data(FILE *file) {
    free_data(file, nullptr);
}

void free_data(unsigned char *pix_data) {
    free_data(nullptr, pix_data);
}

namespace color {
    void write_header(FILE *file_out, int width, int height, unsigned int max_value) {
        fseek(file_out, 0, SEEK_SET);
        fprintf(file_out, "P%c\n%i %i\n%i\n", '6', width, height, max_value);
    }

    void write_to_file(FILE *file_out, int width, int height, unsigned int max_value,
                       unsigned char *pix_data) {
        color::write_header(file_out, width, height, max_value);
        write_data(file_out, 3 * width * height, pix_data);
    }
}

namespace gray {
    void write_header(FILE *file_out, int width, int height, unsigned int max_value) {
        fseek(file_out, 0, SEEK_SET);
        fprintf(file_out, "P%c\n%i %i\n%i\n", '5', width, height, max_value);
    }

    void write_to_file(FILE *file_out, int width, int height, unsigned int max_value,
                       unsigned char *pix_data) {
        gray::write_header(file_out, width, height, max_value);
        write_data(file_out, width * height, pix_data);
    }
}

void write_to_file(FILE *file_out, char char_header, int width, int height, unsigned int max_value,
                   unsigned char *pix_data) {
    write_header(file_out, char_header, width, height, max_value);
    write_data(file_out, width * height, pix_data);
}

// brightness in [0..255] scale
unsigned char limit_brightness(double brightness) {
    return (unsigned char) std::min(255.0, std::max(0.0, brightness));
}

double from_sRGB(double _brightness) {
    if (_brightness <= 0.0031308) {
        return 323.0 * _brightness / 25.0;
    } else {
        return (211 * pow(_brightness, 5.0 / 12.0) - 11) / 200.0;
    }
}

double from_sRGB(int brightness) {
    double _brightness = brightness / 255.0;
    return from_sRGB(_brightness);
}

double to_sRGB(double _brightness) {
    if (_brightness <= 0.04045) {
        return 25.0 * _brightness / 323;
    } else {
        return pow((200 * _brightness + 11) / 211.0, 12.0 / 5.0);
    }
}

double to_sRGB(int brightness) {
    double _brightness = brightness / 255.0;
    return to_sRGB(_brightness);
}

double change_pix_gamma_to_print(double _brightness, double gamma) {
    // Gamma correction:
    if (gamma > 0) {
        _brightness = std::pow(_brightness, 1.0 / gamma);
    } else {
        _brightness = from_sRGB(_brightness);
    }

    return _brightness;
}

double change_pix_gamma_to_print(unsigned char pix_data, double gamma) {
    double _brightness = pix_data / 255.0;
    return change_pix_gamma_to_print(_brightness, gamma);
}

double change_pix_gamma_from_file(double _brightness, double gamma) {
    if (gamma > 0) {
        _brightness = std::pow(_brightness, gamma);
    } else {
        _brightness = to_sRGB(_brightness);
    }

    return _brightness;
}

double change_pix_gamma_from_file(unsigned char pix_data, double gamma) {
    double _brightness = pix_data / 255.0;
    return change_pix_gamma_from_file(_brightness, gamma);
}

void draw_pix(unsigned char *pix_data, int width, int x, int y, int brightness, double gamma) {
    double _brightness = brightness / 255.0;


    change_pix_gamma_to_print(_brightness, gamma);

    pix_data[width * y + x] = limit_brightness(255 * _brightness);
}

void draw_pix(unsigned char *pix_data, int width, int x, int y, double _brightness, double gamma) {
    change_pix_gamma_to_print(_brightness, gamma);
    pix_data[width * y + x] = limit_brightness(255 * _brightness);
}

void decode_gamma_from_file(unsigned char *pix_data, int k, double gamma) {
    for (int i = 0; i < k; ++i) {
        *(pix_data + i) = limit_brightness(255 * change_pix_gamma_from_file(*(pix_data + i), gamma));
    }
}
