#ifndef FG_JSMATH_H
#define FG_JSMATH_H

/** \file jsMath.h
 * Helpers to evaluate user expressions with the JavaScript engine of Qt, used where the
 * Python tools evaluated Python code.
 */

#include <QString>

class QJSEngine;
class QJSValue;

namespace fg
    {
/** Make the members of the JavaScript Math object global (sin, sqrt, abs, PI...), and define
 * the Python names pi, e, fabs, True and False, so that expressions written for Python's math
 * module work unchanged. */
void installMathGlobals(QJSEngine &engine);

/** Translate the Python boolean operators and, or, not of expr to &&, ||, !. Other Python
 * constructs are not translated: chained comparisons such as 0 < x < 1 must be written
 * 0 < x && x < 1. */
QString pythonToJs(const QString &expr);

/** If value is a JavaScript error, return its message with the line number, otherwise an
 * empty string. */
QString jsErrorMessage(const QJSValue &value);
    }  // namespace fg

#endif
