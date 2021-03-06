#ifndef SAILINGFORM_H
#define SAILINGFORM_H

#include <QWidget>
#include <QVector2D>
#include <QtMath>
#include <QFileDialog>
#include <QTextStream>
#include <random>
#include <QGraphicsView>

#include "messagedialog.h"
#include "contourrecognition.h"

namespace Ui {
class SailingForm;
}

class SailingForm : public QWidget
{
    Q_OBJECT

public:
    explicit SailingForm(QWidget *parent = 0);
    ~SailingForm();

private slots:
    void on_startPushButton_clicked();

private:
    Ui::SailingForm *ui;
    QVector2D getUnitStepVector(double Nerror, double speed, double navIntervalMin);
    double getNorthError(int time, int okta, int num);
    QMap<int, double> getTimeElevationMap(QString filename);
    int getGaussianRandomNumber(double mu, double sigma, QString mode);
    int getUniformRandomNumber(int low, int high);
    void drawUnitVectors(QImage &image, QColor &color, QList<QVector2D> &vectorList, QPointF shift); //shift x: horizontal, y: vertical
    void drawNavigationEndPoint(QImage &image, QColor &color, QList<QVector2D> &vectorList, QPointF shift);
    void fitImage(QImage &image, QGraphicsView *view);
    double getNavigationIntervalError(int interval);

    void selectVikingRoute(QString inpath, QString outpath);

    double distance;
    QGraphicsScene scene1, scene2;
    double sumLength;
    ContourRecognition contour;
};

#endif // SAILINGFORM_H
