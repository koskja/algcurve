#pragma once
#include "core.hpp"
#include "thread.hpp"
#include <array>
#include <cassert>
#include <fstream>
#include <numeric>
#include <unordered_map>
#include <vector>

/// A base class for pixels. The values are all in the range [0.0, 1.0].
struct Pixel {
    virtual ~Pixel() = default;
    virtual double r() const = 0;
    virtual double g() const = 0;
    virtual double b() const = 0;
    virtual void set_r(double r) = 0;
    virtual void set_g(double g) = 0;
    virtual void set_b(double b) = 0;
};

struct Triplet : Pixel {
    double _r, _g, _b;
    double r() const override {
        return _r;
    }
    double g() const override {
        return _g;
    }
    double b() const override {
        return _b;
    }
    void set_r(double r) override {
        _r = r;
    }
    void set_g(double g) override {
        _g = g;
    }
    void set_b(double b) override {
        _b = b;
    }
};

struct Real : Pixel {
    double _v;
    double r() const override {
        return _v;
    }
    double g() const override {
        return _v;
    }
    double b() const override {
        return _v;
    }
    void set_r(double r) override {
        _v = r;
    }
    void set_g(double g) override {
        _v = g;
    }
    void set_b(double b) override {
        _v = b;
    }
    void set_v(double v) {
        _v = v;
    }
};

struct BlackWhite : Pixel {
    u8 _v = 0;
    double r() const override {
        return _v;
    }
    double g() const override {
        return _v;
    }
    double b() const override {
        return _v;
    }
    void set_r(double r) override {
        _v = r == 0.0 ? 0 : 1;
    }
    void set_g(double g) override {
        _v = g == 0.0 ? 0 : 1;
    }
    void set_b(double b) override {
        _v = b == 0.0 ? 0 : 1;
    }
    void set_v(u8 v) {
        _v = v;
    }
};

template <usize D, typename P = Triplet> struct Texture {
    union {
        std::array<usize, D> dims;
        struct {
            usize width, height, depth;
        };
    };
    std::vector<P> data;
    Texture() = default;
    template <typename... Args> Texture(Args... args) {
        static_assert(sizeof...(Args) == D, "Number of arguments must match texture dimensions");
        dims = {static_cast<usize>(args)...};
        data.resize(size());
    }

    usize size() const {
        return std::accumulate(dims.begin(), dims.end(), (usize)1, std::multiplies<usize>());
    }

    template <typename... Args> constexpr bool in_bounds(Args... _args) const {
        std::array<usize, D> args = {static_cast<usize>(_args)...};
        for (usize i = 0; i < D; ++i) {
            if (args[i] >= dims[i]) {
                return false;
            }
        }
        return true;
    }

    template <typename... Args> constexpr usize to_index(Args... _args) const {
        std::array<usize, D> args = {static_cast<usize>(_args)...};
        usize index = 0;
        usize multiplier = 1;
        for (usize i = 0; i < D; ++i) {
            index += args[i] * multiplier;
            multiplier *= dims[i];
        }
        return index;
    }

    template <typename Iter> void overwrite(Iter begin, Iter end) {
        if (std::distance(begin, end) != size()) {
            throw std::runtime_error("Iterator range size does not match image dimensions");
        }
        std::copy(begin, end, data.begin());
    }

    Texture<D + 1, P> expand(usize ndimsize, std::function<void(P, std::span<P>)> func) {
        Texture<D + 1, P> result;
        result.dims[0] = ndimsize;
        for (usize i = 0; i < D; ++i) {
            result.dims[i + 1] = dims[i];
        }
        result.data.resize(result.size());
        parallel_for(
            result.size(), [&](usize i) { func(result.data[i], result.data.subspan(i * ndimsize, ndimsize)); }, 1);
        return result;
    }

    template <typename... Args> P& operator()(Args... args) {
        assert(in_bounds(args...));
        return data[to_index(args...)];
    }

    template <typename... Args> const P& operator()(Args... args) const {
        assert(in_bounds(args...));
        return data[to_index(args...)];
    }

    static constexpr u32 info_header_size = 40;
    static constexpr u32 file_header_size = 14;

    template <typename T> void write(std::ostream& os, T value) {
        os.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

    void write_info_header(std::ostream& os) {
        write<u32>(os, info_header_size); // sizeof BITMAPINFOHEADER
        write<u32>(os, (u32)width);
        write<u32>(os, (u32)height);
        write<u16>(os, 1);   // planes
        write<u16>(os, 24);  // bits per pixel
        write<u32>(os, 0);   // compression
        write<u32>(os, 0);   // image size; 0 for uncompressed
        write<u32>(os, 512); // x resolution
        write<u32>(os, 512); // y resolution
        write<u32>(os, 0);   // colors used; default
        write<u32>(os, 0);   // important colors; default
    }

    void write_bmp(std::ostream& os) {
        os.write("BM", 2);
        write<u32>(os, file_header_size + info_header_size + data.size() * (usize)3); // file size
        write<u32>(os, 0);                                                     // reserved
        write<u32>(os, file_header_size + info_header_size);                   // raster data offset
        write_info_header(os);
        for (isize y = height - 1; y >= 0; --y) {
            for (isize x = 0; x < width; ++x) {
                auto& pixel = (*this)(x, y);
                write<u8>(os, pixel.b() * 255);
                write<u8>(os, pixel.g() * 255);
                write<u8>(os, pixel.r() * 255);
            }
        }
    }

    void save_bmp(std::string_view filename) {
        std::ofstream os(filename.data(), std::ios::binary);
        if (!os) {
            throw std::runtime_error("Failed to open file for writing: " + std::string(filename));
        }
        write_bmp(os);
        if (!os) {
            throw std::runtime_error("Failed to write BMP data");
        }
        os.close();
        if (!os) {
            throw std::runtime_error("Failed to close file");
        }
    }
};

template <usize D, typename P = Triplet> struct SparseTexture {
    std::array<usize, D> dims;
    std::unordered_map<usize, P> data;

    SparseTexture() = default;
    template <typename... Args> SparseTexture(Args... args) {
        static_assert(sizeof...(Args) == D, "Number of arguments must match texture dimensions");
        dims = {static_cast<usize>(args)...};
    }

    template <typename... Args> constexpr usize to_index(Args... _args) const {
        std::array<usize, D> args = {static_cast<usize>(_args)...};
        usize index = 0;
        usize multiplier = 1;
        for (usize i = 0; i < D; ++i) {
            index += args[i] * multiplier;
            multiplier *= dims[i];
        }
        return index;
    }

    template <typename... Args> P& operator()(Args... args) {
        return data[to_index(args...)];
    }

    Texture<D, P> to_dense() const {
        Texture<D, P> dense_texture;
        dense_texture.dims = dims;
        dense_texture.data.resize(dense_texture.size());

        for (const auto& [index, pixel] : data) {
            dense_texture.data[index] = pixel;
        }
        return dense_texture;
    }
};

template <typename P = Triplet> using Texture1D = Texture<1, P>;

template <typename P = Triplet> using Texture2D = Texture<2, P>;

template <typename P = Triplet> using Texture3D = Texture<3, P>;
