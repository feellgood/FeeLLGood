#include "jsMath.h"

#include <QJSEngine>
#include <QJSValue>
#include <QRegularExpression>

namespace fg
    {
void installMathGlobals(QJSEngine &engine)
    {
    engine.evaluate(
            "(function (g) {"
            "  Object.getOwnPropertyNames(Math).forEach(function (k) { g[k] = Math[k]; });"
            "  g.pi = Math.PI; g.e = Math.E; g.fabs = Math.abs;"
            "  g.True = true; g.False = false;"
            "})(this);");
    }

QString pythonToJs(const QString &expr)
    {
    QString s = expr;
    s.replace(QRegularExpression("\\band\\b"), "&&");
    s.replace(QRegularExpression("\\bor\\b"), "||");
    s.replace(QRegularExpression("\\bnot\\b"), "!");
    return s;
    }

QString jsErrorMessage(const QJSValue &value)
    {
    if (!value.isError())
        {
        return QString();
        }
    return QString("line %1: %2").arg(value.property("lineNumber").toInt()).arg(value.toString());
    }
    }  // namespace fg
