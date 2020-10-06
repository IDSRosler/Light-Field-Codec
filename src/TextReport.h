//
// Created by cmbos on 9/11/2020.
//

#ifndef LF_CODEC_TEXTREPORT_H
#define LF_CODEC_TEXTREPORT_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <sstream>

struct TextReport {

    explicit TextReport();
    explicit TextReport(std::ostream& os);

    void header(const std::initializer_list<std::string>& columns);
    template <typename T>
    void set_key(const std::string&& key, const T& value);
    void set_width(int width);
    void set_separator(const std::string& sep);
    void print();
private:
    std::unordered_map<std::string, std::string> m_row;
    std::vector<std::string> m_columns;
    std::ostream &m_os;
    int m_width;
    std::string m_sep;
};

TextReport::TextReport() : TextReport(std::cout) { }
TextReport::TextReport(std::ostream &os) : m_os(os), m_width(0), m_sep("|") { }

void TextReport::header(const std::initializer_list<std::string>& columns) {
    m_columns = columns;
    for (const auto& col: m_columns)
        m_width = std::max(m_width, static_cast<int>(col.size() + 3));
    //m_os << m_sep;
    for (const auto& col: m_columns) {
        m_os << std::setw(m_width);
        m_os << col;
        //  m_os << m_sep;
    }
    m_os << std::endl << std::endl;
    m_os.flush();
}
template <typename T>
void TextReport::set_key(const std::string&& key, const T& value) {
    std::stringstream ss;
    ss << value;
    m_row[key] = ss.str();
}
template<>
void TextReport::set_key(const std::string&& key, const std::string& value) {
    m_row[key] = value;
}
void TextReport::print() {
    m_os << m_sep;
    static const std::string placeholder = "-";
    for (const auto& col: m_columns) {
        // TODO: Make m_width customizable for each column.
        m_os << std::setw(m_width);
        auto el = m_row.find(col);
        m_os << (el != m_row.end() ? (el->second) : placeholder);
        m_os << m_sep;
    }
    m_os << std::endl;
    m_row.clear();
    m_os.flush();
}


#endif //LF_CODEC_TEXTREPORT_H
