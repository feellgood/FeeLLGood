#include "quiver2DWidget.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include <QFontMetrics>
#include <QPainter>
#include <QPainterPath>
#include <QPolygonF>

namespace
    {
const QColor ARROW_COLOR(0, 0, 255);
}

std::vector<double> niceTicks(double lo, double hi)
    {
    std::vector<double> ticks;
    const double range = hi - lo;
    if (!(range > 0) || !std::isfinite(range))
        {
        ticks.push_back(lo);
        return ticks;
        }
    const double rough = range / 5.0;
    const double mag = std::pow(10.0, std::floor(std::log10(rough)));
    const double f = rough / mag;
    const double step = (f < 1.5 ? 1.0 : f < 3.5 ? 2.0 : f < 7.5 ? 5.0 : 10.0) * mag;
    for (double t = std::ceil(lo / step) * step; t <= hi + 1e-9 * range; t += step)
        {
        ticks.push_back(std::abs(t) < 1e-12 * range ? 0.0 : t);
        }
    return ticks;
    }

Quiver2DWidget::Quiver2DWidget(std::vector<double> x_, std::vector<double> y_,
                               std::vector<double> u_, std::vector<double> v_,
                               const QString &title_, QWidget *parent)
    : QWidget(parent), x(std::move(x_)), y(std::move(y_)), u(std::move(u_)), v(std::move(v_)),
      title(title_)
    {
    const auto [xlo, xhi] = std::minmax_element(x.begin(), x.end());
    const auto [ylo, yhi] = std::minmax_element(y.begin(), y.end());
    double big = std::max(*yhi - *ylo, *xhi - *xlo);
    if (!(big > 0))
        {
        big = 1.0;
        }
    const double xc = 0.5 * (*xlo + *xhi), yc = 0.5 * (*ylo + *yhi);
    xMin = xc - big;
    xMax = xc + big;
    yMin = yc - big;
    yMax = yc + big;
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);
    }

void Quiver2DWidget::paintEvent(QPaintEvent *)
    {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    p.fillRect(rect(), Qt::white);
    const QFontMetrics fm(font());

    // square plot area, leaving room for the labels and the title
    const int left = fm.horizontalAdvance("-0.000000") + 12, bottom = fm.height() + 12;
    const int top = 2 * fm.height() + 8, right = 20;
    const int side = std::max(10, std::min(width() - left - right, height() - top - bottom));
    const QRectF area(left + (width() - left - right - side) / 2.0,
                      top + (height() - top - bottom - side) / 2.0, side, side);

    auto toX = [&](double xv) { return area.left() + (xv - xMin) / (xMax - xMin) * side; };
    auto toY = [&](double yv) { return area.bottom() - (yv - yMin) / (yMax - yMin) * side; };

    // frame, ticks and labels
    p.setPen(QPen(Qt::black, 1));
    p.drawRect(area);
    for (double t : niceTicks(xMin, xMax))
        {
        const double px = toX(t);
        p.drawLine(QPointF(px, area.bottom()), QPointF(px, area.bottom() + 5));
        const QString s = QString::number(t, 'g', 6);
        p.drawText(QPointF(px - fm.horizontalAdvance(s) / 2.0, area.bottom() + 6 + fm.ascent()),
                   s);
        }
    for (double t : niceTicks(yMin, yMax))
        {
        const double py = toY(t);
        p.drawLine(QPointF(area.left() - 5, py), QPointF(area.left(), py));
        const QString s = QString::number(t, 'g', 6);
        p.drawText(QPointF(area.left() - 8 - fm.horizontalAdvance(s), py + fm.ascent() / 2.0 - 1),
                   s);
        }
    QFont titleFont = font();
    titleFont.setPointSizeF(titleFont.pointSizeF() * 1.2);
    p.setFont(titleFont);
    p.drawText(QRectF(area.left(), 0, side, top), Qt::AlignCenter, title);

    // arrows, scaled as matplotlib does by default
    const std::size_t n = x.size();
    if (n == 0)
        {
        return;
        }
    double sumLength = 0;
    for (std::size_t i = 0; i < n; ++i)
        {
        sumLength += std::hypot(u[i], v[i]);
        }
    const double meanLength = sumLength / n;
    if (!(meanLength > 0))
        {
        return;
        }
    const double sn = std::max(10.0, std::sqrt(static_cast<double>(n)));
    const double pixelsPerUnit = side / (1.8 * meanLength * sn);
    const double sw = std::max(8.0, std::min(25.0, std::sqrt(static_cast<double>(n))));
    const double shaftWidth = 0.06 * side / sw;
    const double headWidth = 3 * shaftWidth, headLength = 5 * shaftWidth;
    const double headAxisLength = 4.5 * shaftWidth;

    p.setClipRect(area);
    p.setPen(Qt::NoPen);
    p.setBrush(ARROW_COLOR);
    for (std::size_t i = 0; i < n; ++i)
        {
        const double length = std::hypot(u[i], v[i]) * pixelsPerUnit;
        if (length <= 0)
            {
            continue;
            }
        // arrow along +x, from 0 to len, shrunk as a whole if shorter than its head
        const double k = std::min(1.0, length / headLength);
        const double hw = headWidth * k, hl = headLength * k, ha = headAxisLength * k;
        const double w = shaftWidth * k;
        const QPolygonF arrow({{0, -w / 2},
                               {length - ha, -w / 2},
                               {length - hl, -hw / 2},
                               {length, 0},
                               {length - hl, hw / 2},
                               {length - ha, w / 2},
                               {0, w / 2}});
        p.save();
        p.translate(toX(x[i]), toY(y[i]));
        p.rotate(-std::atan2(v[i], u[i]) * 180.0 / M_PI);
        p.translate(-length / 2, 0);  // pivot in the middle
        p.drawPolygon(arrow);
        p.restore();
        }
    }
