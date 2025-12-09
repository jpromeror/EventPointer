# Script para comparar archivos modificados uno por uno

$changedFiles = @(
    "ArraysData.R",
    "CDFfromGTF.R",
    "CDFfromGTF_Multipath.R",
    "CreateExSmatrix.R",
    "EventDetection.R",
    "EventDetectionMultipath.R",
    "EventDetection_transcriptome.R",
    "EventPointer.R",
    "EventPointer_Bootstraps.R",
    "EventPointer_IGV.R",
    "EventPointer_RNASeq.R",
    "EventPointer_RNASeq_IGV.R",
    "EventPointer_RNASeq_TranRef.R",
    "EventPointer_RNASeq_TranRef_IGV.R",
    "Events_ReClassification.R",
    "EventXtrans.R",
    "FindPrimers.R",
    "Fit.R",
    "getbootstrapdata.R",
    "GetPSI_FromTranRef.R",
    "Protein_Domain_Enrichment.R",
    "PSI_Statistic.R",
    "ResulTable.R",
    "SF_Prediction.R"
)

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Quedan $($changedFiles.Count) archivos por revisar" -ForegroundColor Yellow
Write-Host "========================================`n" -ForegroundColor Cyan

$index = 1
foreach ($file in $changedFiles) {
    Write-Host "[$index/$($changedFiles.Count)] $file" -ForegroundColor Green
    Write-Host "  IZQUIERDA: EP_Cesar/R/$file (nuevo)" -ForegroundColor Cyan
    Write-Host "  DERECHA: R/$file (actual)`n" -ForegroundColor Yellow
    
    # Abrir en modo comparación
    code -d "EP_Cesar\R\$file" "R\$file"
    
    Write-Host "Presiona Enter cuando hayas terminado con este archivo (o 'q' para salir): " -NoNewline
    $continue = Read-Host
    
    if ($continue -eq 'q') {
        Write-Host "`nSaliendo... Archivos restantes: $($changedFiles.Count - $index)" -ForegroundColor Red
        return
    }
    
    Write-Host ""
    $index++
}

Write-Host "========================================" -ForegroundColor Green
Write-Host "¡Todos los archivos revisados!" -ForegroundColor Green
Write-Host "========================================`n" -ForegroundColor Green
