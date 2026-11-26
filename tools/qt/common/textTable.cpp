#include "textTable.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace fg
    {
TextTable readTextTable(const std::string &fileName)
    {
    std::ifstream f(fileName);
    if (!f)
        {
        throw std::runtime_error("cannot open file '" + fileName + "'");
        }
    TextTable table;
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(f, line))
        {
        ++lineNumber;
        if (!line.empty() && line.back() == '\r')
            {
            line.pop_back();
            }
        const std::size_t first = line.find_first_not_of(" \t");
        if (first == std::string::npos)
            {
            continue;
            }
        if (line[first] == '#')
            {
            table.comments.push_back(line);
            continue;
            }
        std::istringstream iss(line);
        std::vector<double> row;
        std::string field;
        while (iss >> field)
            {
            char *end = nullptr;
            const double x = std::strtod(field.c_str(), &end);
            if (end == field.c_str() || *end != '\0')
                {
                throw std::runtime_error(fileName + ":" + std::to_string(lineNumber)
                                         + ": invalid number '" + field + "'");
                }
            row.push_back(x);
            }
        table.rows.push_back(std::move(row));
        }
    return table;
    }

bool metadataValue(const std::string &line, const std::string &key, std::string &value)
    {
    if (line.compare(0, key.size(), key) != 0)
        {
        return false;
        }
    std::size_t start = key.size();
    if (start < line.size() && (line[start] == ' ' || line[start] == '\t'))
        {
        ++start;
        }
    value = line.substr(start);
    return true;
    }
    }  // namespace fg
