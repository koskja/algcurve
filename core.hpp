#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <span>

typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef size_t usize;
typedef ptrdiff_t isize;
typedef float f32;
typedef double f64;

typedef i16 exp_t;

#define SIMD_SIZE (512 / 8)
#define SIMD_ALIGN SIMD_SIZE
#define SIMD_NUM_VALUES (SIMD_SIZE / sizeof(double))

#ifdef _WIN32
#define __aligned_alloc(align, size) _aligned_malloc(size, align)
#define __aligned_free(ptr) _aligned_free(ptr)
#else
#define __aligned_alloc(align, size) std::aligned_alloc(align, size)
#define __aligned_free(ptr) std::free(ptr)
#endif

template <typename T, usize ALIGN> struct SimdHeapArray {
    T *data;
    usize byte_size;
    usize public_size;
    SimdHeapArray() : data(nullptr), byte_size(0), public_size(0) {}
    SimdHeapArray(usize size) {
        auto bytes = size * sizeof(T);
        auto aligned_bytes = (bytes + (ALIGN - 1)) & ~(ALIGN - 1);
        byte_size = aligned_bytes;
        public_size = size;
        if (aligned_bytes == 0) {
            data = nullptr;
        } else {
            data = (T *)__aligned_alloc(ALIGN, aligned_bytes);
            for (usize i = 0; i < size; ++i) {
                new (&data[i]) T{};
            }
        }
    }
    usize size() const {
        return public_size;
    }
    SimdHeapArray(const SimdHeapArray& other) {
        if (other.byte_size == 0) {
            data = nullptr;
            byte_size = 0;
            public_size = 0;
        } else {
            byte_size = other.byte_size;
            public_size = other.public_size;
            data = (T *)__aligned_alloc(ALIGN, other.byte_size);
            for (usize i = 0; i < other.size(); ++i) {
                new (&data[i]) T{other.data[i]};
            }
        }
    }
    void destroy_data() {
        if (!data) return;
        for (usize i = 0; i < public_size; ++i) {
            data[i].~T();
        }
        __aligned_free(data);
        data = nullptr;
        byte_size = 0;
        public_size = 0;
    }
    SimdHeapArray(SimdHeapArray&& other) noexcept
        : data(other.data), byte_size(other.byte_size), public_size(other.public_size) {
        other.data = nullptr;
        other.byte_size = 0;
        other.public_size = 0;
    }
    SimdHeapArray& operator=(const SimdHeapArray& other) {
        if (this == &other) return *this;
        destroy_data();
        if (other.byte_size == 0) {
            return *this;
        }
        byte_size = other.byte_size;
        public_size = other.public_size;
        data = (T *)__aligned_alloc(ALIGN, other.byte_size);
        for (usize i = 0; i < other.size(); ++i) {
            new (&data[i]) T{other.data[i]};
        }
        return *this;
    }
    SimdHeapArray& operator=(SimdHeapArray&& other) {
        if (this == &other) {
            return *this;
        }
        destroy_data();
        data = other.data;
        byte_size = other.byte_size;
        public_size = other.public_size;
        other.data = nullptr;
        other.byte_size = 0;
        other.public_size = 0;
        return *this;
    }
    ~SimdHeapArray() {
        destroy_data();
    }
    T& operator[](usize idx) {
        return data[idx];
    }

    const T& operator[](usize idx) const {
        return data[idx];
    }

    std::span<T> slice_len(usize index, usize len) {
        return std::span<T>(data + index, len);
    }
    std::span<const T> slice_len(usize index, usize len) const {
        return std::span<const T>(data + index, len);
    }

    std::span<T> slice(usize start, usize end) {
        return std::span<T>(data + start, data + end);
    }
    std::span<const T> slice(usize start, usize end) const {
        return std::span<const T>(data + start, data + end);
    }
};