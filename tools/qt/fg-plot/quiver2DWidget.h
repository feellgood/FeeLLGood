#ifndef FGPLOT_QUIVER2DWIDGET_H
#define FGPLOT_QUIVER2DWIDGET_H

/** \file quiver2DWidget.h
 * 2D quiver plot drawn with QPainter.
 */

#include <QString>
#include <QWidget>
#include <vector>

/** \class Quiver2DWidget
 * Square 2D plot of arrows (u, v) centered on the points (x, y), with axes, ticks and a title.
 * The arrows are scaled as matplotlib's quiver does by default: the mean arrow length is
 * 1 / (1.8 max(10, sqrt(N))) of the plot width, N being the number of arrows.
 */
class Quiver2DWidget : public QWidget
    {
public:
    /** Build the plot of the arrows (u[i], v[i]) at (x[i], y[i]); the four vectors must have
     * the same length. The plotted area is centered on the data and its half width is the
     * largest extent of the data along x or y. */
    Quiver2DWidget(std::vector<double> x, std::vector<double> y, std::vector<double> u,
                   std::vector<double> v, const QString &title, QWidget *parent = nullptr);

    /** Preferred size of the window. */
    QSize sizeHint() const override { return {640, 640}; }

protected:
    /** Draw the frame, the ticks with their labels, the title and the arrows. */
    void paintEvent(QPaintEvent *event) override;

private:
    std::vector<double> x, y, u, v;
    QString title;
    double xMin, xMax, yMin, yMax;  //!< limits of the plotted area
    };

/** Tick positions in [lo, hi] at a "nice" step (1, 2 or 5 times a power of ten) giving about
 * five ticks. */
std::vector<double> niceTicks(double lo, double hi);

#endif
