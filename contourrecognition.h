#ifndef CONTOURRECOGNITION_H
#define CONTOURRECOGNITION_H

#include <QObject>
#include <QWidget>
#include <QtMath>
#include <QVector2D>
#include <QVector3D>
#include <QMap>
#include <QPainter>
#include <QFile>
#include <QTextStream>

class ContourRecognition : public QObject
{
    Q_OBJECT
public:
    explicit ContourRecognition(QObject *parent = 0);
    ~ContourRecognition();

    QList<QPointF> blowUpContour(QList<QPointF> &dataList, double blow, QImage &image);
    void skeletonize(int posX, int posY, QImage &image);
    QList<QPointF> scaleDataToImage(QString dataPath, QImage &image);
    QList<QPointF> scaleContour(QString dataPath, QImage &image);
    double blowDistance(double R, double s, double H, double h);
    QList<QPointF>getRelativeContourPositions(QImage &image);
    QPointF rescaleDataToMapPoints(QPointF data);

private:
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0, height = 0, width = 0;
};

#endif // CONTOURRECOGNITION_H
