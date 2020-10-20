#ifndef __UTILS_H__
#define __UTILS_H__

#include "Typedef.h"
#include "Point4D.h"


#include <vector>
#include <string>
#include <sstream>

#include <opencv2/opencv.hpp>
#include <cmath>
#include <iomanip>
#include <fstream>




#define FOREACH_4D_IDX(idx, shape, stride)                                               \
    for (std::remove_const<decltype(shape.v)>::type v = 0; v < shape.v; v++)             \
        for (std::remove_const<decltype(shape.u)>::type u = 0; u < shape.u; u++)         \
            for (std::remove_const<decltype(shape.y)>::type y = 0; y < shape.y; y++)     \
                for (std::remove_const<decltype(shape.x)>::type x = 0; x < shape.x; x++) \
                    for (int _once = 1, idx = calc_offset(x, y, u, v, stride); _once; _once--)
/*
#define FOREACH_4D(ptr, start_ptr, shape, stride)                                 \
    for (int v = 0; v < shape.v; v++)                                             \
        for (int u = 0; u < shape.u; u++)                                         \
            for (int y = 0; y < shape.y; y++)                                     \
                for (int x = 0; x < shape.x; x++)                                 \
                    for (int _once = 1, _index = calc_offset(x, y, u, v, stride); \
                         ({ptr = start_ptr + _index;_once });                     \
                         _once--)

*/
struct PartitionDescriptor
{
    Point4D offset;
    Point4D shape;
    int transform;

    friend std::ostream &operator<<(std::ostream &os, const PartitionDescriptor &d)
    {
        os << "{offset=" << d.offset
           << ", shape=" << d.shape
           << ", transform=" << d.transform << "}";
        return os;
    }
};

struct quadtree
{
    using node_type = std::shared_ptr<quadtree>;
    node_type up_left;
    node_type up_right;
    node_type down_left;
    node_type down_right;
    int transform;
    char partition;

    static node_type make_node() {
        return std::make_shared<quadtree>();
    }

    std::string repr() const
    {
        std::stringstream ss;
        if (up_left == nullptr && up_right == nullptr && down_left == nullptr && down_right == nullptr)
        {
            ss << transform;
        }
        else
        {
            ss << partition;
            ss << up_left->repr();
            ss << up_right->repr();
            ss << down_left->repr();
            ss << down_right->repr();
        }
        std::string str = ss.str();
        return str;
    }



    void set_child(std::size_t child, const node_type& ptr)
    {
        switch (child) {
            case 0:
                up_left = ptr;
                break;
            case 1:
                up_right = ptr;
                break;
            case 2:
                down_left = ptr;
                break;
            case 3:
                down_right = ptr;
                break;

        }
    }
};


template <typename T>
inline Point4D make_stride(T shape);
template <typename size_type>
auto make_shapes(const std::vector<size_type> &from_shape,
                 const std::vector<size_type> &base_shape);
template <typename size_type>
auto make_shapes(const std::vector<size_type> &from_shape, int levels);

inline int calc_offset(Point4D::value_type x, Point4D::value_type y, Point4D::value_type u, Point4D::value_type v, const Point4D &stride);
inline int calc_offset(const Point4D &point, const Point4D &stride);
template <typename T, typename size_type>
void segment_block(const T *from_block,
                   const std::vector<size_type> &from_shape,
                   T *into_block,
                   int levels);
template <typename T, typename size_type>
void join_segments(const T *from_block,
                   const std::vector<size_type> &from_shape,
                   T *into_block,
                   int levels);
template <typename T1, typename T2>
inline void display_report(std::ostream &os, const T1 &left, const T2 &right);
template <typename T>
std::string to_human_readable(T size);
template <typename T>
std::vector<index_t> calculate_rank(std::vector<T> array);
template <typename T>
double measure_sortness(std::vector<T> array);
inline void extend_borders(float *block, const Point4D &shape, const Point4D &stride);
void show_block(int channel, float *block, const Point4D &shape, const Point4D &stride, const char *window);
void save_microimage(std::string path, Point4D pos, int channel, float *block, const Point4D &shape, const Point4D &stride, std::string suffix, unsigned flags = 0);
void progress_bar(double progress, int bar_length);
void flip_axis(float *block, unsigned to_flip, unsigned flat_size, Point4D shape, Point4D stride);
std::vector<index_t> generate_scan_order(const Point4D &shape, Point4D &stride);
std::vector<index_t> generate_z_order_curve(const Point4D &shape, Point4D &stride);
std::vector<PartitionDescriptor> parse_descriptor(std::string_view descriptor, const Point4D &_shape, bool unsafe = false);
// std::vector<std::shared_ptr<tree_node>> generate_full_binary_trees(std::size_t total_nodes);
std::vector<std::string> convert_fbt_to_descriptor(std::string tree_repr, std::size_t index);
std::unique_ptr<std::deque<std::pair<Point4D, Point4D>>> split_coordinate(char type, const Point4D& _offset, const Point4D& _shape, bool unsafe = false);
bool is_valid_descriptor(const std::string& descriptor);

template <typename Iter, typename Func>
void for_each_4d(const Iter *block, const Point4D& shape, const Point4D& stride, Func f) {
    for (Point4D::value_type v = 0; v < shape.v; v++)
        for (Point4D::value_type u = 0; u < shape.u; u++)
            for (Point4D::value_type y = 0; y < shape.y; y++)
                for (Point4D::value_type x = 0; x < shape.x; x++)
                {
                    auto offset = calc_offset(x, y, u, v, stride);
                    f(block[offset]);
                }
}

template <typename T>
inline Point4D make_stride(T shape)
{
    Point4D stride;
    stride.x = 1;
    stride.y = shape[0];
    stride.u = stride.y * shape[1];
    stride.v = stride.u * shape[2];
    return stride;
}

template <typename size_type>
auto make_shapes(const std::vector<size_type> &from_shape,
                 const std::vector<size_type> &base_shape)
{
    std::vector<std::vector<size_type>> shapes;
    std::vector<size_type> shape;
    shape.resize(4);
    for (int v = 0; v < from_shape[3]; v += shape[3])
    {
        shape[3] = std::min(base_shape[3], from_shape[3] - v);
        for (int u = 0; u < from_shape[2]; u += shape[2])
        {
            shape[2] = std::min(base_shape[3], from_shape[2] - u);
            for (int y = 0; y < from_shape[1]; y += shape[1])
            {
                shape[1] = std::min(base_shape[1], from_shape[1] - y);
                for (int x = 0; x < from_shape[0]; x += shape[0])
                {
                    shape[0] = std::min(base_shape[0], from_shape[0] - x);
                    shapes.push_back(shape);
                }
            }
        }
    }
    return shapes;
}

template <typename size_type>
auto make_shapes(const std::vector<size_type> &from_shape, int levels)
{

    std::vector<std::vector<size_type>> shapes;
    std::vector segment_shape = from_shape;
    std::vector remainder = from_shape;

    auto it_seg = segment_shape.begin();
    auto it_rem = remainder.begin();
    for (; it_seg < segment_shape.end(); it_seg++, it_rem++)
    {
        *it_seg = *it_seg >> levels;
        *it_rem = *it_rem - (*it_seg << levels);
    }
    for (int v = 0; v < (1 << levels); v++)
    {
        for (int u = 0; u < (1 << levels); u++)
        {
            for (int y = 0; y < (1 << levels); y++)
            {
                for (int x = 0; x < (1 << levels); x++)
                {
                    std::vector curr_seg = segment_shape;
                    if (x < remainder[0])
                        curr_seg[0]++;
                    if (y < remainder[1])
                        curr_seg[1]++;
                    if (u < remainder[2])
                        curr_seg[2]++;
                    if (v < remainder[3])
                        curr_seg[3]++;
                    shapes.push_back(curr_seg);
                }
            }
        }
    }
    return shapes;
}


inline int calc_offset(Point4D::value_type x, Point4D::value_type y, Point4D::value_type u, Point4D::value_type v, const Point4D &stride)
{
    return x * stride.x + y * stride.y + u * stride.u + v * stride.v;
}

inline int calc_offset(const Point4D &point, const Point4D &stride)
{
    return point.x * stride.x + point.y * stride.y + point.u * stride.u + point.v * stride.v;
}

template <typename T, typename size_type>
void segment_block(const T *from_block,
                   const std::vector<size_type> &from_shape,
                   T *into_block,
                   int levels)
{
    auto *curr_seg = into_block;
    auto stride = make_stride(from_shape);
    int seg_index = 0;
    const auto shapes = make_shapes(from_shape, levels);
    auto shape = shapes[seg_index];
    for (int v = 0; v < from_shape[3]; v += shape[3])
    {
        for (int u = 0; u < from_shape[2]; u += shape[2])
        {
            for (int y = 0; y < from_shape[1]; y += shape[1])
            {
                for (int x = 0; x < from_shape[0]; x += shape[0])
                {
                    shape = shapes[seg_index++];
                    auto into_stride = make_stride(shape);
                    for (int dv = 0; dv < shape[3]; dv++)
                        for (int du = 0; du < shape[2]; du++)
                            for (int dy = 0; dy < shape[1]; dy++)
                            {
                                auto from_offset = stride.x * x + stride.y * (y + dy) +
                                                   stride.u * (u + du) + stride.v * (v + dv);

                                std::copy(from_block + from_offset,
                                          from_block + from_offset + shape[0], curr_seg);
                                curr_seg += shape[0];
                            }
                }
            }
        }
    }
}
template <typename T, typename size_type>
void join_segments(const T *from_block,
                   const std::vector<size_type> &from_shape,
                   T *into_block,
                   int levels)
{
    auto *curr_seg = from_block;
    auto stride = make_stride(from_shape);
    int seg_index = 0;
    const auto shapes = make_shapes(from_shape, levels);
    auto shape = shapes[seg_index];
    for (int v = 0; v < from_shape[3]; v += shape[3])
    {
        for (int u = 0; u < from_shape[2]; u += shape[2])
        {
            for (int y = 0; y < from_shape[1]; y += shape[1])
            {
                for (int x = 0; x < from_shape[0]; x += shape[0])
                {
                    shape = shapes[seg_index++];
                    auto into_stride = make_stride(shape);
                    for (int dv = 0; dv < shape[3]; dv++)
                        for (int du = 0; du < shape[2]; du++)
                            for (int dy = 0; dy < shape[1]; dy++)
                            {
                                auto into_offset = stride.x * x + stride.y * (y + dy) +
                                                   stride.u * (u + du) + stride.v * (v + dv);
                                std::copy(curr_seg, curr_seg + shape[0], into_block + into_offset);
                                curr_seg += shape[0];
                            }
                }
            }
        }
    }
}

template <typename T1, typename T2>
inline void display_report(std::ostream &os, const T1 &left, const T2 &right)
{
    os << std::setfill('.') << std::setw(30) << std::left << left;
    os << " ";
    os << std::setfill(' ') << std::setw(30) << std::left << right;
    os << std::endl;
}

template <typename T>
std::string to_human_readable(T size)
{
    std::stringstream ss;
    if (size > ((1 << 30) - 1))
        ss << ((double)size / ((1 << 30) - 1)) << "GB";
    else if (size > ((1 << 20) - 1))
        ss << ((double)size / ((1 << 20) - 1)) << "MB";
    else if (size > ((1 << 10) - 1))
        ss << ((double)size / ((1 << 10) - 1)) << "KB";
    else
        ss << size << "B";
    return ss.str();
}



template <typename T>
std::vector<index_t> calculate_rank(std::vector<T> array)
{
    std::map<T, int> rank_map;
    std::vector<T> sorted(array.size());
    std::vector<index_t> ranks;
    int rank = 1;

    std::partial_sort_copy(array.begin(), array.end(), sorted.begin(), sorted.end(),
                           [](auto &a, auto &b) -> bool { return std::abs(a) > std::abs(b); });
    for (auto &e : sorted)
    {
        auto search = rank_map.find(e);
        if (search == rank_map.end())
            rank_map[e] = rank++;
    }

    ranks.reserve(array.size());
    for (auto &e : array)
    {
        auto search = rank_map.find(e);
        ranks.push_back(search->second);
    }

    return ranks;
}

template <typename T>
double measure_sortness(std::vector<T> array)
{
    double array_mean = 0;
    double ranks_mean = 0;
    double a = 0, b = 0, c = 0;

    auto rank_vector = calculate_rank(array);
    for (int i = 0; i < array.size(); i++)
    {
        array_mean += array[i];
        ranks_mean += rank_vector[i];
    }
    array_mean /= (double)array.size();
    ranks_mean /= (double)array.size();

    for (int i = 0; i < array.size(); i++)
    {
        a += ((double)array[i] - array_mean) * ((double)rank_vector[i] - ranks_mean);
        b += ((double)array[i] - array_mean) * ((double)array[i] - array_mean);
        c += ((double)rank_vector[i] - ranks_mean) * ((double)rank_vector[i] - ranks_mean);
    }
    return a / (std::pow(b, 2) * std::pow(c, 2));
}

std::tuple<Point4D, Point4D> partition_at(const Point4D &from_shape,
                  const Point4D &from_offset,
                  size_t at_position,
                  char partition_type);
#endif // __UTILS_H__