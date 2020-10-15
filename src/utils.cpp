#include "utils.h"
#include "Transform.h"
#include <algorithm> // std::sort
#include <chrono>
#include <ctime>
#include <unordered_map>
#include <utility>
#include <string>
#include <stack>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <memory>

void show_block(int channel, float *block, const Point4D &shape,
                const Point4D &stride, const char *window) {

    const char *channel_name;
    cv::String title = "Block Viewer ";
    auto size = cv::Size(13 * 15, 13 * 15);
    cv::Mat mat = cv::Mat::zeros(size, CV_32F);

    for (Point4D::value_type v = 0; v < shape.v; v++)
        for (Point4D::value_type u = 0; u < shape.u; u++)
            for (Point4D::value_type y = 0; y < shape.y; y++)
                for (Point4D::value_type x = 0; x < shape.x; x++) {
                    int index = x * stride.x + y * stride.y + u * stride.u + v * stride.v;
                    auto mat_x = x * shape.u + u;
                    auto mat_y = y * shape.v + v;
                    // mat.at<float>(mat_y, mat_x) = (block[index] + 512.0F) / 1024.0F;
                    auto abs_value = std::abs(block[index]);
                    auto bitsize = std::ceil(std::log2(abs_value));
                    mat.at<float>(mat_y, mat_x) = bitsize / 20;
                }

    switch (channel) {
        case 0:
            channel_name = "Y";
            break;
        case 1:
            channel_name = "Cb";
            break;
        case 2:
            channel_name = "Cr";
            break;
    }
    title += channel_name;
    title += " ";
    title += window;
    cv::namedWindow(title, cv::WINDOW_AUTOSIZE);
    cv::imshow(title, mat);
    cv::waitKey(0);

}

void save_microimage(std::string path, Point4D pos, int channel, float *block, const Point4D &shape,
                     const Point4D &stride, std::string suffix) {
    auto size = cv::Size(shape.x * shape.u, shape.y * shape.v);
    cv::Mat mat = cv::Mat::zeros(size, CV_32F);

    FOREACH_4D_IDX(index, shape, stride) {
                        auto mat_x = x * shape.u + u;
                        auto mat_y = y * shape.v + v;
                        // mat.at<float>(mat_y, mat_x) = 255 * (block[index] + 512.0F) / 1024.0F;
                        auto abs_value = std::abs(block[index]);
                        auto bitsize = std::ceil(std::log2(abs_value));
                        mat.at<float>(mat_y, mat_x) = (bitsize / 20.0) * 255;
                    }

    std::stringstream ss;

    ss << path
       << "v" << pos.v / 13
       << "_u" << pos.u / 13
       << "_y" << pos.y / 15
       << "_x" << pos.x / 15;
    switch (channel) {
        case 0:
            ss << "_Y_";
            break;
        case 1:
            ss << "_Cb_";
            break;
        case 2:
            ss << "_Cr_";
            break;
    }
    ss << suffix << ".png";
    imwrite(ss.str(), mat);
}

void progress_bar(double progress, int bar_length) {
    static std::chrono::time_point<std::chrono::steady_clock> start;
    if (int(progress * 1000) == 0)
        start = std::chrono::steady_clock::now();

    int64_t eta = 0;
    int64_t elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - start)
            .count();

    if (progress > 0)
        eta = ((1 - progress) * elapsed / progress);
    if (progress == 1) {
        for (int i = 0; i < 80; i++)
            std::cout << " ";
        std::cout << "\r";
    } else {
        std::cout << "[";
        int pos = bar_length * progress;
        for (int i = 0; i < bar_length; ++i) {
            if (i < pos)
                std::cout << "=";
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " % "
                  << "ETA " << eta / 60 << "min " << eta % 60 << "sec \r";
    }

    std::cout.flush();
}

void flip_axis(float *block, unsigned to_flip, unsigned flat_size,
               Point4D shape, Point4D stride) {
    float _block[flat_size];
    for (Point4D::value_type v = 0; v < shape.v; v++) {
        for (Point4D::value_type u = 0; u < shape.u; u++) {
            for (Point4D::value_type y = 0; y < shape.y; y++) {
                for (Point4D::value_type x = 0; x < shape.x; x++) {
                    auto dx = to_flip & Transform::AXIS_X ? shape.x - 1 - 2 * x : 0;
                    auto dy = to_flip & Transform::AXIS_Y ? shape.y - 1 - 2 * y : 0;
                    auto du = to_flip & Transform::AXIS_U ? shape.u - 1 - 2 * u : 0;
                    auto dv = to_flip & Transform::AXIS_V ? shape.v - 1 - 2 * v : 0;
                    auto f_offset = calc_offset(x, y, u, v, stride);
                    auto r_offset = calc_offset(x + dx, y + dy, u + du, v + dv, stride);
                    _block[r_offset] = block[f_offset];
                }
            }
        }
    }
    for (Point4D::value_type v = 0; v < shape.v; v++) {
        for (Point4D::value_type u = 0; u < shape.u; u++) {
            for (Point4D::value_type y = 0; y < shape.y; y++) {
                for (Point4D::value_type x = 0; x < shape.x; x++) {
                    auto index = calc_offset(x, y, u, v, stride);
                    block[index] = _block[index];
                }
            }
        }
    }
}

inline void extend_borders(float *block, const Point4D &shape,
                           const Point4D &stride) {
    for (Point4D::value_type y = 0; y < shape.y; y++)
        for (Point4D::value_type x = 0; x < shape.x; x++)
            for (Point4D::value_type v = 0; v < shape.v; v++) {
                Point4D::value_type u = 0;
                int value;
                // Left extension
                while ((value = block[calc_offset(x, y, u, v, stride)]) == 0)
                    u++;
                for (Point4D::value_type i = 0; i < u; i++)
                    block[calc_offset(x, y, i, v, stride)] = value;

                // Right extension
                u = shape.u - 1;
                while ((value = block[calc_offset(x, y, u, v, stride)]) == 0)
                    u--;
                for (Point4D::value_type i = shape.u - 1; i > u; i--)
                    block[calc_offset(x, y, i, v, stride)] = value;
            }
}

struct PairPoint4DHasher {
    std::size_t operator()(const std::pair<Point4D, Point4D> &pair) const {
        auto[p, q] = pair;
        return p.x ^ ((p.y << 1) ^ ((p.u << 2) ^ (p.v << 3))) ^ q.x ^
               ((q.y << 1) ^ ((q.u << 2) ^ (q.v << 3)));
    }
};

static std::unordered_map<std::pair<Point4D, Point4D>, std::vector<index_t>,
        PairPoint4DHasher>
        scan_order_map;

std::vector<index_t> generate_scan_order(const Point4D &shape, Point4D &stride) {
    auto SIZE = shape.getNSamples();
    auto key = std::pair(shape, stride);
    try {
        return scan_order_map.at(key);
    } catch (...) {
        std::vector<index_t> indexes;
        indexes.reserve(SIZE);

        Point4D p(0, 0, 0, 0);
        while (p != shape) {
            if (p < shape) {
                auto _offset = (index_t) calc_offset(p, stride);
                indexes.push_back(_offset);
            }
            if (p[0] > 0) {
                p[1]++;
                p[0]--;
            } else {
                int i;
                for (i = 1; i < 3; i++) {
                    if (p[i] > 0) {
                        p[i + 1]++;
                        p[0] = p[i] - 1;
                        p[i] = 0;
                        break;
                    }
                }
                if (i == 3) {
                    p[0] = p[3] + 1;
                    p[3] = 0;
                }
            }
        }
        scan_order_map[key] = indexes;
        return indexes;
    }
}

std::vector<index_t> generate_z_order_curve(const Point4D &shape, Point4D &stride) {
    auto SIZE = shape.getNSamples();
    std::vector<index_t> scan_order(SIZE);
    auto it = scan_order.begin();

    FOREACH_4D_IDX(_, shape, stride) {
                        index_t index = 0;
                        decltype(x) shadow[4] = {x, y, u, v};
                        int bit = 0;
                        while (shadow[0] | shadow[1] | shadow[2] | shadow[3]) {
                            for (int i = 0; i < 4; i++) {
                                index |= (shadow[i] & 1) << bit;
                                shadow[i] >>= 1;
                                bit++;
                            }
                        }
                        *it++ = index;
                    }
    return scan_order;
}


std::unique_ptr<std::deque<std::pair<Point4D, Point4D>>>
split_coordinate(char type, const Point4D &_offset, const Point4D &_shape, bool unsafe) {
    std::unique_ptr<std::deque<std::pair<Point4D, Point4D>>> stack = std::make_unique<std::deque<std::pair<Point4D, Point4D>>>();

    Point4D shape = _shape;
    Point4D offset = _offset;
    Point4D new_shape = shape;
    Point4D new_offset = offset;

    std::size_t minimum_size = type == 's'
                               ? EncoderParameters::parameters.transform_min_spatial_size
                               : EncoderParameters::parameters.transform_min_angular_size;
    std::size_t i = (type == 's') ? 0 : 2;

    new_shape[i] >>= 1;
    new_shape[i + 1] >>= 1;

    if (unsafe || (new_shape[i] >= minimum_size && new_shape[i + 1] >= minimum_size)) {
        shape[i] -= new_shape[i];
        shape[i + 1] -= new_shape[i + 1];
        new_offset[i] += shape[i];
        new_offset[i + 1] += shape[i + 1];
        // First, push down-right node
        // Then, push down-left and up-right nodes
        new_shape.updateNSamples();
        new_offset.updateNSamples();
        stack->push_front(std::make_pair(new_shape, new_offset));

        auto old = new_shape[i];
        new_shape[i] = shape[i];
        new_offset[i] -= shape[i];
        new_shape.updateNSamples();
        new_offset.updateNSamples();
        stack->push_front(std::make_pair(new_shape, new_offset));

        new_shape[i] = old;
        new_offset[i] += shape[i];
        new_shape[i + 1] = shape[i + 1];
        new_offset[i + 1] -= shape[i + 1];
        new_shape.updateNSamples();
        new_offset.updateNSamples();
        stack->push_front(std::make_pair(new_shape, new_offset));
        shape.updateNSamples();
        offset.updateNSamples();
        stack->push_front(std::make_pair(shape, offset));
    } else {
        stack = nullptr;
    }
    return stack;
}

std::vector<PartitionDescriptor> parse_descriptor(std::string_view descriptor, const Point4D &_shape, bool unsafe) {
    Point4D shape = _shape;
    Point4D offset(0, 0, 0, 0);
    std::deque<std::pair<Point4D, Point4D>> stack;
    std::vector<PartitionDescriptor> parsed_descriptors;
    std::size_t i;

    for (const auto &c : descriptor) {
        switch (c) {
            case 's':
            case 'a': {
                auto temp_stack = split_coordinate(c, offset, shape, unsafe);
                if (temp_stack != nullptr) {
                    auto front = temp_stack->front();
                    temp_stack->pop_front();
                    shape = front.first;
                    offset = front.second;
                    std::copy(std::begin(*temp_stack), std::end(*temp_stack), std::back_inserter(stack));
                }

            }
                break;
            case '1':
            case '2':
            case '3': {
                PartitionDescriptor desc;
                desc.offset = offset;
                desc.shape = shape;
                desc.transform = static_cast<int>(c - '0');
                parsed_descriptors.push_back(desc);

                if (!stack.empty()) {
                    auto front = stack.front();
                    stack.pop_front();
                    shape = front.first;
                    offset = front.second;
                }
            }
        }
    }
    return parsed_descriptors;
}

//
//std::vector<std::shared_ptr<tree_node>> generate_full_binary_trees(std::size_t total_nodes) {
//    std::vector<std::shared_ptr<tree_node>> tree;
//    if (total_nodes == 1) {
//        tree.push_back(std::make_shared<tree_node>());
//    } else if ((total_nodes - 1) % 4 == 0) {
//        for (std::size_t x1 = 1; x1 < total_nodes - 1; x1 += 4)
//            for (std::size_t x2 = 1; x2 < total_nodes - x1 - 1; x2 += 4)
//                for (std::size_t x3 = 1; x3 < total_nodes - x2 - x1 - 1; x3 += 4) {
//                    auto ul_nodes = generate_full_binary_trees(x1); // up left
//                    auto ur_nodes = generate_full_binary_trees(x2); // up right
//                    auto dl_nodes = generate_full_binary_trees(x3); // down left
//                    auto dr_nodes = generate_full_binary_trees(total_nodes - x3 - x2 - x1 - 1); // down right
//                    for (const auto &ul : ul_nodes)
//                        for (const auto &ur : ur_nodes)
//                            for (const auto &dl : dl_nodes)
//                                for (const auto &dr:  dr_nodes) {
//                                    auto root = std::make_shared<tree_node>();
//                                    root->up_left = ul;
//                                    root->up_right = ur;
//                                    root->down_left = dl;
//                                    root->down_right = dr;
//                                    tree.push_back(root);
//                                }
//                }
//    }
//    return tree;
//}
//
//std::vector<std::string> convert_fbt_to_descriptor(std::string tree_repr, std::size_t index) {
//    std::vector<std::string> converted;
//    std::vector<std::string> prev;
//    std::string T_CHOICES;
//    std::transform(std::cbegin(Transform::T_CHOICES),
//                   std::cend(Transform::T_CHOICES),
//                   std::back_inserter(T_CHOICES),
//                   [](auto choice) {
//                       return '0' + static_cast<unsigned char>(choice);
//                   });
//
//    if (index == tree_repr.size()) {
//        converted.push_back(tree_repr);
//    } else {
//        std::string prefix(tree_repr.substr(0, index > 0 ? index : 0));
//        std::string suffix(tree_repr.substr(index + 1, tree_repr.size() - 1));
//        std::string choices;
//
//        switch (tree_repr[index]) {
//            case 'P':
//                choices = Transform::P_CHOICES;
//                break;
//            case 'T':
//                choices = T_CHOICES;
//                break;
//            default:
//                choices = "";
//                prev = convert_fbt_to_descriptor(tree_repr, index + 1);
//                std::copy(std::begin(prev), std::end(prev), std::back_inserter(converted));
//        }
//
//        for (const auto &c : choices) {
//            std::stringstream ss;
//            ss << prefix << c << suffix;
//            prev = convert_fbt_to_descriptor(ss.str(), index + 1);
//            std::copy(std::begin(prev), std::end(prev), std::back_inserter(converted));
//        }
//    }
//    return converted;
//}
//


//def verify_descriptor(desc):
//    count = 1
//    for c in desc:
//        if count == 0:
//            return False
//
//        if c in {'P', 's', 'a'}:
//            count += 3
//        elif c in {'T', '1', '2', '3'}:
//            count -= 1
//        else:
//            return False
//
//    else:
//        return count == 0

bool is_valid_descriptor(const std::string &descriptor) {
    int count = 1;
    for (const auto &c: descriptor) {
        if (count == 0)
            return false;

        switch (c) {
            case 'P':
            case 's':
            case 'a':
                count += 3;
                break;

            case 'T':
            case '1':
            case '2':
            case '3':
                count--;
                break;

            default:
                return false;
        }
    }
    return count == 0;
}



std::tuple<Point4D, Point4D> partition_at(const Point4D &from_shape,
                                          const Point4D &from_offset,
                                          size_t at_position,
                                          char partition_type)
{
    Point4D shape = from_shape;
    Point4D offset = from_offset;
    size_t s, t; // Axis for partition

    switch (partition_type)
    {
        case 'S':
        case 's':
            s = 0;
            t = 1;
            break;
        case 'A':
        case 'a':
            s = 2;
            t = 3;
            break;
        default:
            throw std::invalid_argument("partition_type expects either 's' or 'a'. ");
    }

    switch (at_position)
    {
        case 0:
            shape[s] = (shape[s] + 1) >> 1U;
            shape[t] = (shape[t] + 1) >> 1U;
            break;
        case 1:
            shape[s] = shape[s] >> 1U;
            shape[t] = (shape[t] + 1) >> 1U;
            offset[s] += (from_shape[s] + 1) >> 1U;
            break;
        case 2:
            shape[s] = (shape[s] + 1) >> 1U;
            shape[t] = shape[t] >> 1U;
            offset[t] += (from_shape[t] + 1) >> 1U;
            break;
        case 3:
            shape[s] = shape[s] >> 1U;
            shape[t] = shape[t] >> 1U;
            offset[s] += (from_shape[s] + 1) >> 1U;
            offset[t] += (from_shape[t] + 1) >> 1U;
            break;
        default:
            throw std::invalid_argument("at_position expects values between 0 and 3, inclusive. ");
    }
    return std::make_tuple(shape, offset);
}
