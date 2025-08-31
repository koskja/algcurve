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

template <typename T, usize ALIGN> struct SimdHeapArray {
    T *data;
    usize byte_size;
    SimdHeapArray() : data(nullptr), byte_size(0) {}
    SimdHeapArray(usize size) {
        auto bytes = size * sizeof(T);
        auto aligned_bytes = (bytes + (ALIGN - 1)) & ~(ALIGN - 1);
        byte_size = aligned_bytes;
        if (aligned_bytes == 0) {
            data = nullptr;
        } else {
            data = (T *)aligned_alloc(ALIGN, aligned_bytes);
        }
    }
    SimdHeapArray(const SimdHeapArray& other) {
        if (other.byte_size == 0) {
            data = nullptr;
            byte_size = 0;
        } else {
            byte_size = other.byte_size;
            data = (T *)aligned_alloc(ALIGN, other.byte_size);
            memcpy(data, other.data, other.byte_size);
        }
    }
    SimdHeapArray(SimdHeapArray&& other) noexcept : data(other.data), byte_size(other.byte_size) {
        other.data = nullptr;
        other.byte_size = 0;
    }
    SimdHeapArray& operator=(const SimdHeapArray& other) {
        if (this == &other) return *this;
        if (data) free(data);
        if (other.byte_size == 0) {
            data = nullptr;
            byte_size = 0;
        } else {
            byte_size = other.byte_size;
            data = (T *)aligned_alloc(ALIGN, other.byte_size);
            memcpy(data, other.data, other.byte_size);
        }
        return *this;
    }
    SimdHeapArray& operator=(SimdHeapArray&& other) {
        if (this == &other) {
            return *this;
        }
        if (data) free(data);
        data = other.data;
        byte_size = other.byte_size;
        other.data = nullptr;
        other.byte_size = 0;
        return *this;
    }
    ~SimdHeapArray() {
        if (data) free(data);
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
        return std::span<T>(data + start, end - start);
    }
    std::span<const T> slice(usize start, usize end) const {
        return std::span<const T>(data + start, end - start);
    }
};