#ifndef CONTOURRECOGNITION_H
#define CONTOURRECOGNITION_H

#include <QObject>
#include <QWidget>

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

signals:

public slots:
};

#endif // CONTOURRECOGNITION_H
