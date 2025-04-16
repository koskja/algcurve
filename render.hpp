#include "polynomial.hpp"
#include "image.hpp"

#define SIMD_SIZE 512 / 8
#define SIMD_ALIGN SIMD_SIZE
#define SIMD_NUM_VALUES (SIMD_SIZE / sizeof(double))

struct ImageParams {
    std::function<std::array<double, 2>(usize, usize)> to_plane;
    std::function<double(double)> soft_clamp;
    usize width, height;
};

struct AnimationParams {
    ImageParams image;
    Polynomial<double, 2> p1, p2;
    std::span<double> lambdas;
};

template <typename DYNAMIC, typename STATIC, typename LEN=u64, typename INDEX=LEN>
class CompressedArray {
    std::vector<INDEX> m_indices;
    std::vector<STATIC> m_static_data;
    std::vector<DYNAMIC> m_dynamic_data;
    public:
    std::span<INDEX> indices() {
        return std::span(m_indices);
    }
    std::span<STATIC> static_data() {
        return std::span(m_static_data);
    }
    std::span<DYNAMIC> dynamic_data() {
        return std::span(m_dynamic_data);
    }
    usize size() {
        return m_indices.size();
    }
    LEN length(usize i) {
        assert(i < m_indices.size());
        if (i == m_indices.size() - 1) return m_dynamic_data.size() - m_indices[i];
        return m_indices[i + 1] - m_indices[i];
    }
    const std::span<DYNAMIC> operator[](usize i) const {
        assert(i < m_indices.size());
        usize length = length(i);
        return std::span(m_dynamic_data.data() + m_indices[i], length);
    }
    std::span<DYNAMIC> operator[](usize i) {
        assert(i < m_indices.size());
        usize length = length(i);
        return std::span(m_dynamic_data.data() + m_indices[i], length);
    }
    const STATIC& operator()(usize i) const {
        assert(i < m_static_data.size());
        return m_static_data[i];
    }
    STATIC& operator()(usize i) {
        assert(i < m_static_data.size());
        return m_static_data[i];
    }
    void push_back(DYNAMIC data) {
        m_dynamic_data.push_back(data);
    }
    void push_back(STATIC data) {
        usize true_size = m_dynamic_data.size();
        auto index = static_cast<INDEX>(true_size);
        assert(index == true_size && "Array has more elements than INDEX can cover");
        usize true_len = m_dynamic_data.size() - m_indices.back();
        assert(true_len == (LEN)true_len && "Subarray has more elements than LEN can cover");
        m_static_data.push_back(data);
        m_indices.push_back(m_dynamic_data.size());
    }
    template <typename T>
    void push_back(STATIC data, T dynamic_data) {
        auto span = std::span<DYNAMIC>(dynamic_data);
        m_static_data.push_back(data);
        m_indices.push_back(m_dynamic_data.size());
        m_dynamic_data.insert(m_dynamic_data.end(), dynamic_data.begin(), dynamic_data.end());
    }
    template<typename U, typename F>
    std::vector<U> flatten_map(F _func) {
        std::vector<U> result;
        std::function<U(STATIC&, DYNAMIC&)> func = _func;
        result.resize(m_dynamic_data.size());
        parallel_for(m_dynamic_data.size(), [&](usize i) {
            auto index = m_indices[i];
            parallel_for(length(i), [&](usize j) {
                result[index + j] = func(m_static_data[i], m_dynamic_data[index + j]);
            });
        });
        return result;
    }
};