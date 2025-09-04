#pragma once
#include "core.hpp"
#include <array>
#include <cassert>
#include <fstream>
#include <vector>

struct PaletteColor {
    u8 index;

    static constexpr u8 BLACK = 0;
    static constexpr u8 BLUE = 1;
    static constexpr u8 GREEN = 2;
    static constexpr u8 AQUA = 3;
    static constexpr u8 RED = 4;
    static constexpr u8 PURPLE = 5;
    static constexpr u8 YELLOW = 6;
    static constexpr u8 WHITE = 7;
    static constexpr u8 GRAY = 8;
    static constexpr u8 LIGHT_BLUE = 9;
    static constexpr u8 LIGHT_GREEN = 10;
    static constexpr u8 LIGHT_AQUA = 11;
    static constexpr u8 LIGHT_RED = 12;
    static constexpr u8 LIGHT_PURPLE = 13;
    static constexpr u8 LIGHT_YELLOW = 14;
    static constexpr u8 BRIGHT_WHITE = 15;
};

struct RGB {
    u8 r, g, b;
};

constexpr std::array<RGB, 16> a_palette = {{
    {0, 0, 0},       // BLACK
    {0, 0, 128},     // BLUE
    {0, 128, 0},     // GREEN
    {0, 128, 128},   // AQUA
    {128, 0, 0},     // RED
    {128, 0, 128},   // PURPLE
    {128, 128, 0},   // YELLOW
    {192, 192, 192}, // WHITE
    {128, 128, 128}, // GRAY
    {0, 0, 255},     // LIGHT_BLUE
    {0, 255, 0},     // LIGHT_GREEN
    {0, 255, 255},   // LIGHT_AQUA
    {255, 0, 0},     // LIGHT_RED
    {255, 0, 255},   // LIGHT_PURPLE
    {255, 255, 0},   // LIGHT_YELLOW
    {255, 255, 255}  // BRIGHT_WHITE
}};

struct Image {
    Image() : _width(0), _height(0), _data(0) {}
    Image(usize width, usize height) : _width(width), _height(height), _data(width * height, {PaletteColor::BLACK}) {}

    PaletteColor& operator()(usize x, usize y) {
        assert(x < _width && y < _height);
        return _data[y * _width + x];
    }

    const PaletteColor& operator()(usize x, usize y) const {
        assert(x < _width && y < _height);
        return _data[y * _width + x];
    }

    usize width() const {
        return _width;
    }
    usize height() const {
        return _height;
    }

    void save_bmp(std::string_view filename) {
        std::ofstream os(filename.data(), std::ios::binary);
        if (!os) {
            throw std::runtime_error("Failed to open file for writing: " + std::string(filename));
        }

        std::vector<u8> compressed_data = rle8_compress();

        const u32 file_header_size = 14;
        const u32 info_header_size = 40;
        const u32 palette_size = 256 * 4;
        const u32 pixel_data_offset = file_header_size + info_header_size + palette_size;
        const u32 file_size = pixel_data_offset + compressed_data.size();

        // File Header
        os.write("BM", 2);
        write<u32>(os, file_size);
        write<u32>(os, 0); // Reserved
        write<u32>(os, pixel_data_offset);

        // Info Header
        write<u32>(os, info_header_size);
        write<u32>(os, (u32)_width);
        write<u32>(os, (u32)_height);
        write<u16>(os, 1); // Planes
        write<u16>(os, 8); // Bits per pixel
        write<u32>(os, 1); // BI_RLE8 compression
        write<u32>(os, (u32)compressed_data.size());
        write<u32>(os, 2835); // X pixels per meter (72 DPI)
        write<u32>(os, 2835); // Y pixels per meter (72 DPI)
        write<u32>(os, 256);  // Colors used
        write<u32>(os, 256);  // Important colors

        // Color Table
        for (const auto& color : a_palette) {
            write<u8>(os, color.b);
            write<u8>(os, color.g);
            write<u8>(os, color.r);
            write<u8>(os, 0); // Reserved
        }
        // Pad palette to 256 colors
        for (usize i = a_palette.size(); i < 256; ++i) {
            write<u8>(os, 0);
            write<u8>(os, 0);
            write<u8>(os, 0);
            write<u8>(os, 0);
        }

        // Pixel Data
        os.write(reinterpret_cast<const char *>(compressed_data.data()), compressed_data.size());

        if (!os) {
            throw std::runtime_error("Failed to write BMP data");
        }
    }

  private:
    usize _width, _height;
    std::vector<PaletteColor> _data;

    template <typename T> void write(std::ostream& os, T value) {
        os.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

    std::vector<u8> rle8_compress() const {
        std::vector<u8> compressed;
        for (isize y = _height - 1; y >= 0; --y) {
            std::vector<u8> row_data;
            isize x = 0;
            while (x < _width) {
                isize run_length = 1;
                while (x + run_length < _width && run_length < 255 &&
                       (*this)(x + run_length, y).index == (*this)(x, y).index) {
                    run_length++;
                }

                u8 index = (*this)(x, y).index;
                row_data.push_back(run_length);
                row_data.push_back(index);
                x += run_length;
            }
            row_data.push_back(0); // End of line
            row_data.push_back(0);
            compressed.insert(compressed.end(), row_data.begin(), row_data.end());
        }

        compressed.push_back(0); // End of bitmap
        compressed.push_back(1);

        return compressed;
    }
};
