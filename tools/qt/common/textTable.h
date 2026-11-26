#ifndef FG_TEXTTABLE_H
#define FG_TEXTTABLE_H

/** \file textTable.h
 * Reading of the text files feeLLGood writes, .sol and .evol: comment lines starting with '#'
 * and rows of numbers separated by tabs or spaces.
 */

#include <string>
#include <vector>

namespace fg
    {
/** Content of a .sol or .evol file. */
struct TextTable
    {
    std::vector<std::string> comments;      //!< lines starting with '#', without end of line
    std::vector<std::vector<double>> rows;  //!< numbers of the other non empty lines
    };

/** Read a .sol or .evol file. Empty lines are skipped. Throws std::runtime_error if the file
 * cannot be opened or if a field of a data line is not a number. */
TextTable readTextTable(const std::string &fileName);

/** If line starts with key (for instance "## time:"), store in value what follows the key,
 * without the leading blank, and return true. Otherwise return false. */
bool metadataValue(const std::string &line, const std::string &key, std::string &value);
    }  // namespace fg

#endif
