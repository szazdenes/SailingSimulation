#ifndef CONTOURRECOGNITION_H
#define CONTOURRECOGNITION_H

#include <QObject>
#include <QWidget>
#include <QtMath>
#include <QVector2D>
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

    QImage getContour(QString imagePath);
    QImage blowUpContour(QImage &image);
    QPointF blownPoint(double x, double y, double blow);
    void skeletonize(int posX, int posY, QImage &image);
    QList<QPoint> getNeighbourList(QImage &image);

    QList<QPointF> scaleDataToImage(QString dataPath, QImage &image);
    double blowDistance(double R, double s, double H, double h);

signals:

public slots:
};

#endif // CONTOURRECOGNITION_H
