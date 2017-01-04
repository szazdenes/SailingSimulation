#include "sailingform.h"
#include "ui_sailingform.h"

SailingForm::SailingForm(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SailingForm)
{
    ui->setupUi(this);

    ui->trajectoryGraphicsView->setScene(&scene1);
    ui->multipleRunGraphicsView->setScene(&scene2);
}

SailingForm::~SailingForm()
{
    delete ui;
}

QVector2D SailingForm::getUnitStepVector(double Nerror, double speed)
{
    QVector2D result;
    result.setX(speed * qCos(Nerror * M_PI / 180.0));
    result.setY(speed * qSin(Nerror * M_PI / 180.0));
    return result;
}

double SailingForm::getNorthError(int time, int okta, int num)
{
    double northError;
    QMap<int, double> elevationMap, roundedElevationMap;
    QFile file;

    if(ui->solRadioButton->isChecked()){
        if(ui->calciteCheckBox->isChecked() && num==1){
            if(time <= 12)
                file.setFileName("../cal_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("../cal_sol_pm_ave.csv");
        }
        else if(ui->cordieriteCheckBox->isChecked() && num==2){
            if(time <= 12)
                file.setFileName("../cord_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("../cord_sol_pm_ave.csv");
        }
        else if(ui->tourmalineCheckBox->isChecked() && num==3){
            if(time <= 12)
                file.setFileName("../tour_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("../tour_sol_pm_ave.csv");
        }
        else
            return -999;
        elevationMap = getTimeElevationMap("../elevation_Bergen_sol.dat");
    }

    if(ui->equRadioButton->isChecked()){
        if(ui->calciteCheckBox->isChecked() && num==1){
            if(time <= 12)
                file.setFileName("../cal_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("../cal_equ_pm_ave.csv");
        }

        else if(ui->cordieriteCheckBox->isChecked() && num==2){
            if(time <= 12)
                file.setFileName("../cord_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("../cord_equ_pm_ave.csv");
        }

        else if(ui->tourmalineCheckBox->isChecked() && num==3){
            if(time <= 12)
                file.setFileName("../tour_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("../tour_equ_pm_ave.csv");
        }
        else
            return -999;
        elevationMap = getTimeElevationMap("../elevation_Bergen_equ.dat");
    }

    foreach(int key, elevationMap.keys()){
        if(ui->solRadioButton->isChecked() && qRound(elevationMap[key]) >= 50)
            roundedElevationMap[key] = 50.0;
        else if(ui->equRadioButton->isChecked() && qRound(elevationMap[key]) >= 25)
            roundedElevationMap[key] = 25.0;
        else
            roundedElevationMap[key] = qRound(elevationMap[key]);
    }

    QMap<QPair<int, int>, double> NErrorMap; /*first:elevation, second:okta*/

    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug("Opening error.");
    }
    QTextStream stream(&file);
    stream.readLine();
    while(!stream.atEnd()){
        QString line = stream.readLine();
        QTextStream linestream(&line);
        int elev, cloud;
        double NError;
        linestream >> elev >> cloud >> NError;
        QPair<int, int> keyPair;
        keyPair.first = elev;
        keyPair.second = cloud;
        NErrorMap[keyPair] = NError;
    }

    int roundedElevation;
    if(roundedElevationMap[time] > 0.0 && roundedElevationMap[time] <= 5.0)
        roundedElevation = 5;
    else if(roundedElevationMap[time] > 5 && roundedElevationMap[time] <= 10)
        roundedElevation = 10;
    else if(roundedElevationMap[time] > 10 && roundedElevationMap[time] <= 15)
        roundedElevation = 15;
    else if(roundedElevationMap[time] > 15 && roundedElevationMap[time] <= 20)
        roundedElevation = 20;
    else if(roundedElevationMap[time] > 20 && roundedElevationMap[time] <= 25)
        roundedElevation = 25;
    else if(roundedElevationMap[time] > 25 && roundedElevationMap[time] <= 30)
        roundedElevation = 30;
    else if(roundedElevationMap[time] > 30 && roundedElevationMap[time] <= 35)
        roundedElevation = 35;
    else if(roundedElevationMap[time] > 35 && roundedElevationMap[time] <= 40)
        roundedElevation = 40;
    else if(roundedElevationMap[time] > 40 && roundedElevationMap[time] <= 45)
        roundedElevation = 45;
    else if(roundedElevationMap[time] > 45 && roundedElevationMap[time] <= 50)
        roundedElevation = 50;
    else return -999;

    QPair<int, int> keyPair(roundedElevation, okta);

    northError = NErrorMap[keyPair];

    return northError;
}

QMap<int, double> SailingForm::getTimeElevationMap(QString filename)
{
    QFile file(filename);
    QMap<int, double> timeElevMap;
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug("Opening error.");
    }
    QTextStream stream(&file);
    stream.readLine();
    while(!stream.atEnd()){
        QString line = stream.readLine();
        QTextStream linestream(&line);
        int time;
        double elev;
        linestream >> time >> elev;
        timeElevMap[time] = elev;
    }
    return timeElevMap;
}

int SailingForm::getGaussianRandomNumber(double mu, double sigma) //0,3 works fine
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (double) X2);
    }

    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);

    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    double result = (mu + sigma * (double) X1);
    double roundedResult = qRound(result);
    if(roundedResult <= -8)
        roundedResult = -8;
    if(roundedResult >= 8)
        roundedResult = 8;

    return roundedResult;
}

int SailingForm::getUniformRandomNumber(int low, int high)
{
    return qrand() % ((high + 1) - low) + low;
}

void SailingForm::drawUnitVectors(QImage &image, QGraphicsScene &scene, QColor &color, QList<QVector2D> &vectorList, QPointF startingPoint, QPointF shift)
{
    QPainter painter(&image);
    painter.setPen(color);

    painter.drawLine(QPointF(shift.x() - startingPoint.x(), startingPoint.y() + shift.y()), QPointF(shift.x() - vectorList.at(0).x(), vectorList.at(0).y() + shift.y()));

    QVector2D fromCurrentVector = vectorList.at(0);
    QVector2D toCurrentVector = vectorList.at(0);

    for(int i = 1; i < vectorList.size(); i++){
        toCurrentVector += vectorList.at(i);
        painter.drawLine(QPointF(shift.x() - fromCurrentVector.x(), fromCurrentVector.y() + shift.y()), QPointF(shift.x() - toCurrentVector.x(), toCurrentVector.y() + shift.y()));
        fromCurrentVector += vectorList.at(i);
    }

    QPen pen;
    pen.setColor(Qt::black);
    pen.setWidth(10);
    painter.setPen(pen);
    painter.drawPoint(shift.x() - startingPoint.x(), startingPoint.y() + shift.y());
    pen.setColor(Qt::gray);
    painter.setPen(pen);
    painter.drawPoint(0, shift.y());

    painter.end();
    scene.clear();
    scene.addPixmap(QPixmap::fromImage(image));
}

void SailingForm::drawNavigationEndPoint(QImage &image, QGraphicsScene &scene, QColor &color, QList<QVector2D> &vectorList, QPointF startingPoint, QPointF shift)
{
    QPainter painter(&image);
    QPen pen;
    pen.setColor(color);
    pen.setWidth(3);
    painter.setPen(pen);

    QVector2D toCurrentVector = vectorList.at(0);

    for(int i = 1; i < vectorList.size(); i++)
        toCurrentVector += vectorList.at(i);

    painter.drawPoint(shift.x() - toCurrentVector.x(), toCurrentVector.y() + shift.y());

    pen.setColor(Qt::black);
    pen.setWidth(10);
    painter.setPen(pen);
    painter.drawPoint(shift.x() - startingPoint.x(), startingPoint.y() + shift.y());
    pen.setColor(Qt::gray);
    painter.setPen(pen);
    painter.drawPoint(0, shift.y());

    painter.end();
//    scene.clear();
    scene.addPixmap(QPixmap::fromImage(image));
}

void SailingForm::on_startPushButton_clicked()
{
    QImage endPointImage(ui->multipleRunGraphicsView->width(), ui->multipleRunGraphicsView->height(), QImage::Format_ARGB32_Premultiplied);
    endPointImage.fill(Qt::white);
    QImage trajectoryImage(ui->trajectoryGraphicsView->width(), ui->trajectoryGraphicsView->height(), QImage::Format_ARGB32_Premultiplied);
    trajectoryImage.fill(Qt::white);
    QColor color;

    for (int z = 1; z <= 3; z++){

        if(z==1)
            color = Qt::red;
        if(z==2)
            color = Qt::green;
        if(z==3)
            color = Qt::blue;

        QList<QVector2D> unitStepVectorList;

        for(int i = 0; i < ui->numOfRunsSpinBox->value(); i++){
            unitStepVectorList.clear();
            int firstOkta = getUniformRandomNumber(0,8);
            int currentOkta;
            int currentTime, startingTime, lengthOfDay;

            if(ui->solRadioButton->isChecked()){
                startingTime = 3;
                lengthOfDay = 17;
            }
            if(ui->equRadioButton->isChecked()){
                startingTime = 6;
                lengthOfDay = 11;
            }

            currentOkta = firstOkta;

            for(int i = 0; i < ui->simLengthSpinBox->value(); i++){
                currentTime = startingTime;
                for(int j = 0; j < lengthOfDay; j++){
                    double NError = getNorthError(currentTime, currentOkta, z);
                    if(NError != -999)
                        unitStepVectorList.append(getUnitStepVector(NError, (ui->trajectoryGraphicsView->width()/((double)ui->simLengthSpinBox->value()*17))));

                    currentTime++;
                    currentOkta += getGaussianRandomNumber(0,3);
                    if(currentOkta <= 0)
                        currentOkta = 0;
                    if(currentOkta >= 8)
                        currentOkta = 8;
                }
            }
            if(!unitStepVectorList.isEmpty())
                drawNavigationEndPoint(endPointImage, scene2, color, unitStepVectorList, QPointF(0, 0), QPointF(ui->multipleRunGraphicsView->width(), ui->multipleRunGraphicsView->height()/2.0));
            QApplication::processEvents();
        }
        if(!unitStepVectorList.isEmpty())
            drawUnitVectors(trajectoryImage, scene1, color, unitStepVectorList, QPointF(0, 0), QPointF(ui->trajectoryGraphicsView->width(), ui->trajectoryGraphicsView->height()/2.0));

    }

    MessageDialog messDialog("Simulation ready");
    messDialog.exec();
}
