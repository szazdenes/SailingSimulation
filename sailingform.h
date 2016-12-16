#ifndef SAILINGFORM_H
#define SAILINGFORM_H

#include <QWidget>
#include <QVector2D>
#include <QtMath>

namespace Ui {
class SailingForm;
}

class SailingForm : public QWidget
{
    Q_OBJECT

public:
    explicit SailingForm(QWidget *parent = 0);
    ~SailingForm();

private:
    Ui::SailingForm *ui;
    QVector2D getUnitStepVector(double Nerror);

    double distance;
};

#endif // SAILINGFORM_H
