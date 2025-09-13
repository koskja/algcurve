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

/// A fixed-size container of data of type `T` that is aligned to `SIMD_ALIGN`.
/// `T` must be default constructible, but may be non-trivial to destroy/copy/move.
template <typename T> class SimdHeapArray {
    T *m_data;
    usize m_byte_size;
    usize m_public_size;

    void destroy_data() {
        if (!m_data) return;
        for (usize i = 0; i < m_public_size; ++i) {
            m_data[i].~T();
        }
        __aligned_free(m_data);
        m_data = nullptr;
        m_byte_size = 0;
        m_public_size = 0;
    }

  public:
    SimdHeapArray() : m_data(nullptr), m_byte_size(0), m_public_size(0) {}
    SimdHeapArray(usize size) {
        auto bytes = size * sizeof(T);
        auto aligned_bytes = (bytes + (SIMD_ALIGN - 1)) & ~(SIMD_ALIGN - 1);
        m_byte_size = aligned_bytes;
        m_public_size = size;
        if (aligned_bytes == 0) {
            m_data = nullptr;
        } else {
            m_data = (T *)__aligned_alloc(SIMD_ALIGN, aligned_bytes);
            for (usize i = 0; i < size; ++i) {
                new (&m_data[i]) T{};
            }
        }
    }
    usize size() const {
        return m_public_size;
    }
    SimdHeapArray(const SimdHeapArray& other) {
        if (other.m_byte_size == 0) {
            m_data = nullptr;
            m_byte_size = 0;
            m_public_size = 0;
        } else {
            m_byte_size = other.m_byte_size;
            m_public_size = other.m_public_size;
            m_data = (T *)__aligned_alloc(SIMD_ALIGN, other.m_byte_size);
            for (usize i = 0; i < other.size(); ++i) {
                new (&m_data[i]) T{other.m_data[i]};
            }
        }
    }
    SimdHeapArray(SimdHeapArray&& other) noexcept
        : m_data(other.m_data), m_byte_size(other.m_byte_size), m_public_size(other.m_public_size) {
        other.m_data = nullptr;
        other.m_byte_size = 0;
        other.m_public_size = 0;
    }
    SimdHeapArray& operator=(const SimdHeapArray& other) {
        if (this == &other) return *this;
        destroy_data();
        if (other.m_byte_size == 0) {
            return *this;
        }
        m_byte_size = other.m_byte_size;
        m_public_size = other.m_public_size;
        m_data = (T *)__aligned_alloc(SIMD_ALIGN, other.m_byte_size);
        for (usize i = 0; i < other.size(); ++i) {
            new (&m_data[i]) T{other.m_data[i]};
        }
        return *this;
    }
    SimdHeapArray& operator=(SimdHeapArray&& other) {
        if (this == &other) {
            return *this;
        }
        destroy_data();
        m_data = other.m_data;
        m_byte_size = other.m_byte_size;
        m_public_size = other.m_public_size;
        other.m_data = nullptr;
        other.m_byte_size = 0;
        other.m_public_size = 0;
        return *this;
    }
    ~SimdHeapArray() {
        destroy_data();
    }
    T& operator[](usize idx) {
        return m_data[idx];
    }
    const T& operator[](usize idx) const {
        return m_data[idx];
    }
    std::span<T> slice_len(usize index, usize len) {
        return std::span<T>(m_data + index, len);
    }
    std::span<const T> slice_len(usize index, usize len) const {
        return std::span<const T>(m_data + index, len);
    }
    std::span<T> slice(usize start, usize end) {
        return std::span<T>(m_data + start, m_data + end);
    }
    std::span<const T> slice(usize start, usize end) const {
        return std::span<const T>(m_data + start, m_data + end);
    }
};