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
    QVector2D getUnitStepVector(double Nerror, double speed);
    double getNorthError(int time, int okta, int num);
    QMap<int, double> getTimeElevationMap(QString filename);
    int getGaussianRandomNumber(double mu, double sigma);
    int getUniformRandomNumber(int low, int high);
    void drawUnitVectors(QImage &image, QGraphicsScene &scene, QColor &color, QList<QVector2D> &vectorList, QPointF startingPoint, QPointF shift); //shift x: horizontal, y: vertical
    void drawNavigationEndPoint(QImage &image, QGraphicsScene &scene, QColor &color, QList<QVector2D> &vectorList, QPointF startingPoint, QPointF shift);

    double distance;
    QGraphicsScene scene1, scene2;
};

#endif // SAILINGFORM_H
