#include "sailingform.h"
#include "ui_sailingform.h"

SailingForm::SailingForm(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SailingForm)
{
    ui->setupUi(this);
}

SailingForm::~SailingForm()
{
    delete ui;
}

QVector2D SailingForm::getUnitStepVector(double Nerror)
{
    QVector2D result;
    result.setX(ui->speedDoubleSpinBox->value() * qCos(Nerror * M_PI / 180.0));
    result.setY(ui->speedDoubleSpinBox->value() * qSin(Nerror * M_PI / 180.0));
    return result;
}
