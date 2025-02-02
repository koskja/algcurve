#pragma once
#include <vector>
#include <fstream>
#include "core.hpp"

struct Pixel {
    virtual ~Pixel() = default;
    virtual double r() const = 0;
    virtual double g() const = 0;
    virtual double b() const = 0;
    virtual void set_r(double r) = 0;
    virtual void set_g(double g) = 0;
    virtual void set_b(double b) = 0;
};

struct RGBPixel : Pixel {
    double _r, _g, _b;
    double r() const override { return _r; }
    double g() const override { return _g; }
    double b() const override { return _b; }
    void set_r(double r) override { _r = r; }
    void set_g(double g) override { _g = g; }
    void set_b(double b) override { _b = b; }
};

struct GrayscalePixel : Pixel {
    double _v;
    double r() const override { return _v; }
    double g() const override { return _v; }
    double b() const override { return _v; }
    void set_r(double r) override { _v = r; }
    void set_g(double g) override { _v = g; }
    void set_b(double b) override { _v = b; }
    void set_v(double v) { _v = v; }
};

template<typename P = RGBPixel> requires std::derived_from<P, Pixel>
struct Image {
    usize width;
    usize height;
    std::vector<P> data;

    Image(usize width, usize height) : width(width), height(height), data(width * height) {}
    Image(usize width, usize height, const P& pixel) : width(width), height(height), data(width * height, pixel) {}
    Image(usize width, usize height, std::vector<P> data) : width(width), height(height), data(data) {}
    template<typename Iter>
    Image(usize width, usize height, Iter begin, Iter end) : width(width), height(height) {
        this->overwrite(begin, end);
    }

    template<typename Iter>
    void overwrite(Iter begin, Iter end) {
        if (std::distance(begin, end) != width * height) {
            throw std::runtime_error("Iterator range size does not match image dimensions");
        }
        std::copy(begin, end, data.begin());
    }

    P& operator()(usize x, usize y) {
        return data[y * width + x];
    }

    const P& operator()(usize x, usize y) const {
        return data[y * width + x];
    }

    static constexpr u32 info_header_size = 40;
    static constexpr u32 file_header_size = 14;

    // This is a hack to allow writing to the stream without having to pass it around.
    // It naturally is not thread-safe.
    std::ostream* os = nullptr;
    template<typename T>
    void write(T value) {
        os->write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    void write_info_header() {
        write<u32>(info_header_size); // sizeof BITMAPINFOHEADER
        write<u32>(width);
        write<u32>(height);
        write<u16>(1); // planes
        write<u16>(24); // bits per pixel
        write<u32>(0); // compression
        write<u32>(0); // image size; 0 for uncompressed
        write<u32>(512); // x resolution
        write<u32>(512); // y resolution
        write<u32>(0); // colors used; default
        write<u32>(0); // important colors; default
    }

    void write_bmp() {
        os->write("BM", 2);
        write<u32>(file_header_size + info_header_size + data.size() * 3); // file size
        write<u32>(0); // reserved
        write<u32>(file_header_size + info_header_size); // raster data offset
        write_info_header();
        for (isize y = height - 1; y >= 0; --y) {
            for (isize x = 0; x < width; ++x) {
                auto& pixel = (*this)(x, y);
                write<u8>(pixel.b() * 255);
                write<u8>(pixel.g() * 255);
                write<u8>(pixel.r() * 255);
            }
        }
    }

    void save_bmp(std::string_view filename) {
        std::ofstream os_(filename.data(), std::ios::binary);
        os = &os_;
        if (!os_) {
            throw std::runtime_error("Failed to open file for writing: " + std::string(filename));
        }
        write_bmp();
        if (!os_) {
            throw std::runtime_error("Failed to write BMP data");
        }
        os_.close();
        if (!os_) {
            throw std::runtime_error("Failed to close file");
        }
    }
};


