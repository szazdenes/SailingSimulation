#ifndef SAILINGFORM_H
#define SAILINGFORM_H

#include <QWidget>
#include <QVector2D>
#include <QtMath>
#include <QFileDialog>
#include <QTextStream>
#include <random>
#include <QGraphicsView>

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
    double getNorthError(int time, int okta);
    QMap<int, double> getTimeElevationMap(QString filename);
    int getGaussianRandomNumber(double mu, double sigma);
    int getUniformRandomNumber(int low, int high);
    void drawUnitVectors(QList<QVector2D> &vectorList, QPointF startingPoint, double verticalShift);
    void drawNavigationEndPoint(QList<QVector2D> &vectorList, QPointF startingPoint, double veicalShift);

    double distance;
    QGraphicsScene scene1, scene2;
};

#endif // SAILINGFORM_H
