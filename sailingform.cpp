#include "sailingform.h"
#include "ui_sailingform.h"

SailingForm::SailingForm(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SailingForm)
{
    ui->setupUi(this);

    ui->trajectoryGraphicsView->setScene(&scene1);

    distance = 2720; //km

    ui->tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tableWidget->setHorizontalHeaderLabels(QStringList() << "Good" << "Wrong" << "Success (%)");

}

SailingForm::~SailingForm()
{
    delete ui;
}

QVector2D SailingForm::getUnitStepVector(double Nerror, double speed, double navIntervalMin)
{
    QVector2D result;
    result.setX((speed/60.0)*navIntervalMin * qCos(Nerror * M_PI / 180.0));
    result.setY(-1*(speed/60.0)*navIntervalMin * qSin(Nerror * M_PI / 180.0));
    return result;
}

double SailingForm::getNorthError(int time, int okta, int num)
{
    double northError;
    QMap<int, double> elevationMap, roundedElevationMap;
    QFile file;

    if(ui->solRadioButton->isChecked()){
        if(ui->calciteRadioButton->isChecked() && num==1){
            if((double)time/60.0 < 12)
                file.setFileName("../cal_sol_am.csv");
            if((double)time/60.0 > 12)
                file.setFileName("../cal_sol_pm.csv");
            if((double)time/60.0 == 12){
                int randNum = getUniformRandomNumber(0,1);
                if(randNum == 0)
                    file.setFileName("../cal_sol_am.csv");
                if(randNum == 1)
                    file.setFileName("../cal_sol_pm.csv");
            }
        }
        else if(ui->cordieriteRadioButton->isChecked() && num==2){
            if((double)time/60.0 < 12)
                file.setFileName("../cord_sol_am.csv");
            if((double)time/60.0 > 12)
                file.setFileName("../cord_sol_pm.csv");
            if((double)time/60.0 == 12){
                int randNum = getUniformRandomNumber(0,1);
                if(randNum == 0)
                    file.setFileName("../cord_sol_am.csv");
                if(randNum == 1)
                    file.setFileName("../cord_sol_pm.csv");
            }
        }
        else if(ui->tourmalineRadioButton->isChecked() && num==3){
            if((double)time/60.0 < 12)
                file.setFileName("../tour_sol_am.csv");
            if((double)time/60.0 > 12)
                file.setFileName("../tour_sol_pm.csv");
            if((double)time/60.0 == 12){
                int randNum = getUniformRandomNumber(0,1);
                if(randNum == 0)
                    file.setFileName("../tour_sol_am.csv");
                if(randNum == 1)
                    file.setFileName("../tour_sol_pm.csv");
            }
        }
        else
            return -999;
        elevationMap = getTimeElevationMap("../elevation_Bergen_sol.dat");
    }

    if(ui->equRadioButton->isChecked()){
        if(ui->calciteRadioButton->isChecked() && num==1){
            if((double)time/60.0 < 12)
                file.setFileName("../cal_equ_am.csv");
            if((double)time/60.0 > 12)
                file.setFileName("../cal_equ_pm.csv");
            if((double)time/60.0 == 12){
                int randNum = getUniformRandomNumber(0,1);
                if(randNum == 0)
                    file.setFileName("../cal_equ_am.csv");
                if(randNum == 1)
                    file.setFileName("../cal_equ_pm.csv");
            }
        }

        else if(ui->cordieriteRadioButton->isChecked() && num==2){
            if((double)time/60.0 < 12)
                file.setFileName("../cord_equ_am.csv");
            if((double)time/60.0 > 12)
                file.setFileName("../cord_equ_pm.csv");
            if((double)time/60.0 == 12){
                int randNum = getUniformRandomNumber(0,1);
                if(randNum == 0)
                    file.setFileName("../cord_equ_am.csv");
                if(randNum == 1)
                    file.setFileName("../cord_equ_pm.csv");
            }
        }

        else if(ui->tourmalineRadioButton->isChecked() && num==3){
            if((double)time/60.0 < 12)
                file.setFileName("../tour_equ_am.csv");
            if((double)time/60.0 > 12)
                file.setFileName("../tour_equ_pm.csv");
            if((double)time/60.0 == 12){
                int randNum = getUniformRandomNumber(0,1);
                if(randNum == 0)
                    file.setFileName("../tour_equ_am.csv");
                if(randNum == 1)
                    file.setFileName("../tour_equ_pm.csv");
            }
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

    QMap<QPair<int, int>, QList<double>> NErrorMap; /*first:elevation, second:okta*/

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
        NErrorMap[keyPair].append(NError);
    }

    int roundedElevation;
    if(roundedElevationMap[qRound((double)time/60.0)] >= 0.0 && roundedElevationMap[qRound((double)time/60.0)] <= 5.0)
        roundedElevation = 5;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 5 && roundedElevationMap[qRound((double)time/60.0)] <= 10)
        roundedElevation = 10;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 10 && roundedElevationMap[qRound((double)time/60.0)] <= 15)
        roundedElevation = 15;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 15 && roundedElevationMap[qRound((double)time/60.0)] <= 20)
        roundedElevation = 20;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 20 && roundedElevationMap[qRound((double)time/60.0)] <= 25)
        roundedElevation = 25;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 25 && roundedElevationMap[qRound((double)time/60.0)] <= 30)
        roundedElevation = 30;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 30 && roundedElevationMap[qRound((double)time/60.0)] <= 35)
        roundedElevation = 35;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 35 && roundedElevationMap[qRound((double)time/60.0)] <= 40)
        roundedElevation = 40;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 40 && roundedElevationMap[qRound((double)time/60.0)] <= 45)
        roundedElevation = 45;
    else if(roundedElevationMap[qRound((double)time/60.0)] > 45 && roundedElevationMap[qRound((double)time/60.0)] <= 50)
        roundedElevation = 50;
    else
        return -999;

    QPair<int, int> keyPair(roundedElevation, okta);

    double NErrorNum = getUniformRandomNumber(0, NErrorMap[keyPair].size()-1);
    northError = NErrorMap[keyPair].at(NErrorNum);

    elevation = roundedElevation;
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

int SailingForm::getGaussianRandomNumber(double mu, double sigma, QString mode) //0,3 works fine, mode: nav, cloud
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1)
    {
        call = !call;

        double res = (mu + sigma * (double) X2);
        if(mode == "nav" && res <= 1)
            return 1;
        else
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

    if(mode == "cloud"){
        if(roundedResult <= -8)
            roundedResult = -8;
        if(roundedResult >= 8)
            roundedResult = 8;
    }

    if(mode == "nav"){
        if(roundedResult <= 1)
            roundedResult = 1;
    }

    return roundedResult;
}

int SailingForm::getUniformRandomNumber(int low, int high)
{
    return qrand() % ((high + 1) - low) + low;
}

void SailingForm::drawUnitVectors(QImage &image, QColor &color, QList<QVector2D> &vectorList, QPointF shift)
{
    QPainter painter(&image);
    QPen pen;
    pen.setStyle(Qt::SolidLine);
    pen.setColor(color);
    painter.setPen(pen);

    painter.drawLine(QPointF(shift.x(), shift.y()), QPointF(shift.x() - vectorList.at(0).x(), vectorList.at(0).y() + shift.y()));

    QVector2D fromCurrentVector = vectorList.at(0);
    QVector2D toCurrentVector = vectorList.at(0);

    /*for chechking vector lengths*/
    sumLength = vectorList.at(0).length();
    double sumdifflength = toCurrentVector.length();

    for(int i = 1; i < vectorList.size(); i++){
        toCurrentVector += vectorList.at(i);
        painter.drawLine(QPointF(shift.x() - fromCurrentVector.x(), fromCurrentVector.y() + shift.y()),
                         QPointF(shift.x() - toCurrentVector.x(), toCurrentVector.y() + shift.y()));
        sumdifflength += QVector2D(toCurrentVector-fromCurrentVector).length();
        fromCurrentVector += vectorList.at(i);
        sumLength+=vectorList.at(i).length();
    }

    pen.setColor(Qt::black);
    pen.setWidth(10);
    painter.setPen(pen);
    painter.drawPoint(shift.x(), shift.y());
    pen.setColor(Qt::gray);
    painter.setPen(pen);
    painter.drawPoint((image.width()/86.98)*(-42.7 - (-67.989)), shift.y());

    painter.end();
}

void SailingForm::drawNavigationEndPoint(QImage &image, QColor &color, QList<QVector2D> &vectorList, QPointF shift)
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
    painter.drawPoint(shift.x(), shift.y());
    pen.setColor(Qt::gray);
    painter.setPen(pen);
    painter.drawPoint((image.width()/86.98)*(-42.7 - (-67.989)), shift.y());

    painter.end();
}

void SailingForm::fitImage(QImage &image, QGraphicsView *view)
{

    QRectF imageRect = image.rect();
    QRectF rect = view->viewport()->rect();
    double fitSize = qMin<double>(rect.width() / imageRect.width(), rect.height() / imageRect.height());

    QMatrix matrix = view->matrix().inverted();
    QRectF visibleRect = matrix.mapRect(view->viewport()->rect());
    double zoom = qMin<double>(visibleRect.width() / rect.width(), visibleRect.height() / rect.height());
    zoom *= fitSize;

    view->scale(zoom, zoom);
}

double SailingForm::getNavigationIntervalError(int interval)
{
    double error = (1.0/6.0)*(double)interval;
//    double result = (double)getUniformRandomNumber(0, qRound(error*60))/60.0;
    double result = error;
    int rand = getUniformRandomNumber(0,1);
    if(rand == 0) return result;
    if(rand == 1) return -1*result;
}

void SailingForm::selectVikingRoute(QString inpath, QString outpath)
{
    QFile file(inpath);
    QFile outfile(outpath);

    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
            qDebug("Opening error.");
        }

    if(!outfile.open(QIODevice::WriteOnly | QIODevice::Text)){
            qDebug("Opening error.");
        }
    QTextStream stream(&file), out(&outfile);
    double x, y;

    while(!stream.atEnd()){
        QString line = stream.readLine();
        QTextStream linestream(&line);
        linestream >> x >> y;
        if(x > -68 && x < 19 && y < 84 && y > 50)
            out << x << "\t" << y << "\n";
    }

    file.close();
    outfile.close();
}

void SailingForm::addToList(QString string, bool clear)
{
    if(clear == true) ui->listWidget->clear();
    ui->listWidget->addItem(string);
    ui->listWidget->scrollToBottom();
}

double SailingForm::getMinDistance(QList<QPointF> &contourList, QPointF &currentPoint, QImage &image)
{
    QList<double> distList;
    foreach(QPointF curr, contourList){
        curr.setX(curr.x()*image.width());
        curr.setY(curr.y()*image.height());
        double dist = qSqrt((currentPoint.x()-curr.x())*(currentPoint.x()-curr.x())
                            + (currentPoint.y()-curr.y())*(currentPoint.y()-curr.y()));
        distList.append(dist);
    }
    qSort(distList.begin(), distList.end());
    double minDist = distList.first();
    return minDist;
}

void SailingForm::on_startPushButton_clicked()
{
    QString equsol, stone;
    QString outName;
    QImage background(QSize(800, 600), QImage::Format_ARGB32);
    background.fill(Qt::white);
    QList<QPointF> backgroundPoints = contour.scaleDataToImage("../terkep.dat", background);
    for(int i = 0; i < backgroundPoints.size(); i++){
        background.setPixelColor(QPoint(qRound(backgroundPoints.at(i).x()), qRound(backgroundPoints.at(i).y())), Qt::black);
    }

    QImage trajectoryImage(background.size(), QImage::Format_ARGB32);
    QImage contourImage("../cont_calib.png");
    QImage contAreaImage("../cont_calib_filled.png");
    QImage contImage(background.size(), QImage::Format_ARGB32);
    trajectoryImage.fill(Qt::white);
    contImage.fill(Qt::transparent);
    QColor color;
    QPainter p2(&trajectoryImage);

    double lengthOfVectorList = qAbs((trajectoryImage.width()/86.98)*(5.3 - (-67.989))
                                     - (trajectoryImage.width()/86.98)*(-42.7 - (-67.989)));
    double blowDist = contour.blowDistance(6372797, 1000, 1000, 10); //needs to be scaled
    double blowDistPixel = blowDist*lengthOfVectorList/(distance*1000);

    QList<QPointF> contourPoints = contour.scaleContour("../contour.dat", background);
    //    QList<QPointF> blownContour = contour.blowUpContour(contourPoints, blowDistPixel, background);
    QList<QPointF> relativeContourPoints = contour.getRelativeContourPositions(contourImage);
    QList<QPointF> relativeContAreaPoints = contour.getRelativeContourPositions(contAreaImage);

    p2.drawImage(0, 0, background);
    p2.end();

    fitImage(background, ui->trajectoryGraphicsView);

    double voyageTime = distance/ui->speedDoubleSpinBox->value();

    int startingTime, lengthOfDay;

    for (int z = 1; z <= 3; z++){

        QList<QVector2D> unitStepVectorList;

        int currentStone;
        if(ui->calciteRadioButton->isChecked()){
            stone = "cal";
            currentStone = 1;
        }
        if(ui->cordieriteRadioButton->isChecked()){
            stone = "cord";
            currentStone = 2;
        }
        if(ui->tourmalineRadioButton->isChecked()){
            stone = "tour";
            currentStone = 3;
        }


        if(z==currentStone){

            if(ui->solRadioButton->isChecked()){
                equsol = "sol";
                startingTime = 3;
                lengthOfDay = 17;
            }
            if(ui->equRadioButton->isChecked()){
                equsol = "equ";
                startingTime = 6;
                lengthOfDay = 11;
            }

            QString filename = "trajectory_" + stone + "_" + equsol + "_" + "speed" + QString::number(ui->speedDoubleSpinBox->value()) + "_"
                    + "days" + QString::number(ui->simLengthSpinBox->value()) + "_" + "navperiodicity"
                    + QString::number(ui->hourIntervalSpinBox->value())
                    + "_" + "runs" + QString::number(ui->numOfRunsSpinBox->value()) + ".csv";
            QFile outFile(filename);

            if(!outFile.open(QIODevice::WriteOnly | QIODevice::Text))
                qDebug("baj");
            QTextStream outStream(&outFile);

            //        double resultedNError = -999;
            double good = 0, wrong = 0;
            for(int i = 0; i < ui->numOfRunsSpinBox->value(); i++){
                unitStepVectorList.clear();
                outStream << "----------------------------\n";

                int firstOkta = getUniformRandomNumber(0,8);
                int currentOkta;
                int currentTime;

                currentOkta = firstOkta;

                double NError;
                int navigationInterval;
                double navIntervalWithError;
                int navIntervalWithErrorMin;
                //            double sumNELengthX = 0, sumNELengthY = 0;
                QVector2D endpointVector;
                QPointF shift = QPointF((trajectoryImage.width()/86.98)*(5.3 - (-67.989)),
                                        -1*(trajectoryImage.height()/33.59)*(61 - 83.599));
                bool success = false;
                //            QList<double> minDistList;
                //            QList<double> resNEList;

                for(int k = 0; k < ui->simLengthSpinBox->value(); k++){
                    int counter = 0;
                    currentTime = startingTime*60;
                    navigationInterval = ui->hourIntervalSpinBox->value();
                    navIntervalWithError = navigationInterval;
                    navIntervalWithErrorMin = qRound(60*navIntervalWithError);
                    for(int j = 0; j < lengthOfDay*60; j++){
                        if(counter%navIntervalWithErrorMin == 0){
                            if(ui->noiseCheckBox->isChecked()){
                                NError = getNorthError(currentTime, currentOkta, z);
                                double navError = getNavigationIntervalError(navigationInterval);
                                navIntervalWithError = (double)navigationInterval + navError;
                                navIntervalWithErrorMin = qRound(60*navIntervalWithError);
                                counter = 0;
                            }
                            if(!ui->noiseCheckBox->isChecked()){
                                NError = getNorthError(currentTime, currentOkta, z);
                                counter = 0;
                            }
                            if(NError != -999){
                                //                        sumNELengthX += cos(NError*M_PI/180.0);
                                //                        sumNELengthY += sin(NError*M_PI/180.0);
                                QVector2D unitStepVector = getUnitStepVector(NError, (lengthOfVectorList/voyageTime),
                                                                             navIntervalWithErrorMin);
                                unitStepVectorList.append(unitStepVector); //((double)ui->simLengthSpinBox->value()*17)))); when according sailing days
                                endpointVector += unitStepVector;
                                QPointF unitStepPoint(shift.x() - endpointVector.x(), endpointVector.y() + shift.y());
                                QPointF rescaledPoint = contour.rescaleDataToMapPoints(unitStepPoint);
                                outStream << QString::number(rescaledPoint.x()) << "\t" << QString::number(rescaledPoint.y()) << "\n";
                                double minDist = getMinDistance(relativeContourPoints, unitStepPoint, background);

                                //                        if(sumNELengthY > 0)
                                //                            resultedNError = acos(sumNELengthX/(sqrt(sumNELengthX*sumNELengthX + sumNELengthY*sumNELengthY)))*180.0/M_PI;
                                //                        else
                                //                            resultedNError = -1*acos(sumNELengthX/(sqrt(sumNELengthX*sumNELengthX + sumNELengthY*sumNELengthY)))*180.0/M_PI;

                                //                        resNEList.append(resultedNError);
                                //                        minDistList.append(minDist);

                                if(!relativeContAreaPoints.isEmpty() && minDist < 30 && success == false){
                                    foreach(QPointF currentPixel, relativeContAreaPoints){
                                        if(qRound(currentPixel.x()*background.width()) == qRound(unitStepPoint.x())
                                                && qRound(currentPixel.y()*background.height()) == qRound(unitStepPoint.y())){
                                            success = true;
                                            break;
                                        }
                                    }
                                }

                                //                        qDebug("%f", minDist);


                            }
//                            if(i == 0 && j==0 && NError != -999) addToList("elev: " + QString::number(elevation) + "\t" + "okta: " + QString::number(currentOkta) + "\t" + "NError: " + QString::number(NError), false);
//                            else if(NError != -999) addToList("elev: " + QString::number(elevation) + "\t" + "okta: " + QString::number(currentOkta) + "\t" + "NError: " + QString::number(NError), false);
                        }
                        currentTime++;

                        if(j%60 == 0){
                            currentOkta += getGaussianRandomNumber(0,2, "cloud");
                            if(currentOkta <= 0)
                                currentOkta = 0;
                            if(currentOkta >= 8)
                                currentOkta = 8;
                        }

                        counter++;
                    }
                }

                //            if(!minDistList.isEmpty()){
                //                QList<double> sortedMinDistList = minDistList;
                //                qSort(sortedMinDistList.begin(), sortedMinDistList.end());
                //                double minDistMin = sortedMinDistList.first();
                //                double resNE = resNEList.at(minDistList.indexOf(minDistMin));

                //                if(minDistMin < 5 && resNE > -5 && success == false){
                //                    success = true;
                //                    qDebug("success: %f %f", minDistMin, resNE);
                //                }
                //                qDebug("%f %f", minDistMin, resNE);
                //            }

                if(!success){
                    color = Qt::red;
                    wrong++;
                }
                if(success){
                    color = Qt::green;
                    good++;
                }

                if(!unitStepVectorList.isEmpty()){
                    drawUnitVectors(trajectoryImage, color, unitStepVectorList, shift);
                }
                //            if(resultedNError > -999) addToList("result: " + QString::number(resultedNError), false);
                addToList("run: " + QString::number(i+1), true);
                QApplication::processEvents();
            }


            double successNum = 100*good/(good+wrong);
            ui->tableWidget->setItem(ui->tableWidget->rowCount()-1, 0, new QTableWidgetItem(QString::number(good)));
            ui->tableWidget->setItem(ui->tableWidget->rowCount()-1, 1, new QTableWidgetItem(QString::number(wrong)));
            ui->tableWidget->setItem(ui->tableWidget->rowCount()-1, 2, new QTableWidgetItem(QString::number(successNum)));

            QPainter painter(&trajectoryImage);
            painter.setPen(Qt::magenta);
            for(int i = 0; i < relativeContourPoints.size(); i++){
                contImage.setPixelColor(QPoint(qRound(relativeContourPoints.at(i).x()*background.width()),
                                               qRound(relativeContourPoints.at(i).y()*background.height())), Qt::blue);
                if(i < contourPoints.size()-1){
                    //            background.setPixelColor(QPoint(qRound(blownContour.at(i).x()), qRound(blownContour.at(i).y())), Qt::magenta);
                    painter.drawLine(contourPoints.at(i), contourPoints.at(i+1));
                }
            }
            painter.drawImage(0, 0, contImage);
            painter.end();

            ui->trajectoryGraphicsView->scene()->clear();
            ui->trajectoryGraphicsView->scene()->addPixmap(QPixmap::fromImage(trajectoryImage));

            //    MessageDialog messDialog("Simulation ready");
            //    messDialog.exec();

            outName = "trajectory_" + stone + "_" + equsol + "_" + "speed" + QString::number(ui->speedDoubleSpinBox->value()) + "_"
                    + "days" + QString::number(ui->simLengthSpinBox->value()) + "_" + "navperiodicity"
                    + QString::number(ui->hourIntervalSpinBox->value())
                    + "_" + "runs" + QString::number(ui->numOfRunsSpinBox->value()) + ".png";

            outFile.close();
        }
    }

    trajectoryImage.save(outName);

}
