# Script para comparar archivos entre R/ y EP_Cesar/R/ uno por uno

# Lista de archivos que han cambiado
$changedFiles = @(
    "ArraysData.R",
    "AuxFunctions.R",
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

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Comparador de archivos EventPointer" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Total de archivos modificados: $($changedFiles.Count)" -ForegroundColor Yellow
Write-Host ""

$index = 1
foreach ($file in $changedFiles) {
    Write-Host "[$index/$($changedFiles.Count)] Archivo: $file" -ForegroundColor Green
    Write-Host ""
    Write-Host "Opciones:" -ForegroundColor Yellow
    Write-Host "  1 - Ver diferencias con Compare-Object (líneas diferentes)"
    Write-Host "  2 - Ver diferencias con git diff (formato diff)"
    Write-Host "  3 - Ver estadísticas (líneas añadidas/eliminadas)"
    Write-Host "  4 - Abrir ambos archivos en VS Code"
    Write-Host "  s - Saltar este archivo"
    Write-Host "  q - Salir"
    Write-Host ""
    
    $choice = Read-Host "Elige una opción"
    
    switch ($choice) {
        "1" {
            Write-Host "`n--- Diferencias detectadas ---`n" -ForegroundColor Cyan
            $diff = Compare-Object (Get-Content "R\$file") (Get-Content "EP_Cesar\R\$file") -IncludeEqual:$false
            $diff | Select-Object -First 50 | Format-Table @{Label="Origen";Expression={if($_.SideIndicator -eq "<="){"R/"}else{"EP_Cesar/R/"}}}, InputObject
            if ($diff.Count -gt 50) {
                Write-Host "`n... mostrando solo las primeras 50 diferencias de $($diff.Count) totales`n" -ForegroundColor Yellow
            }
        }
        "2" {
            Write-Host "`n--- Git Diff ---`n" -ForegroundColor Cyan
            git diff --no-index "R\$file" "EP_Cesar\R\$file" | Select-Object -First 100
        }
        "3" {
            $oldLines = (Get-Content "R\$file").Count
            $newLines = (Get-Content "EP_Cesar\R\$file").Count
            $diff = Compare-Object (Get-Content "R\$file") (Get-Content "EP_Cesar\R\$file")
            $removed = ($diff | Where-Object {$_.SideIndicator -eq "<="}).Count
            $added = ($diff | Where-Object {$_.SideIndicator -eq "=>"}).Count
            
            Write-Host "`n--- Estadísticas ---" -ForegroundColor Cyan
            Write-Host "Archivo original (R/): $oldLines líneas" -ForegroundColor White
            Write-Host "Archivo nuevo (EP_Cesar/R/): $newLines líneas" -ForegroundColor White
            Write-Host "Diferencia: $(if($newLines -gt $oldLines){'+'}else{''})$($newLines - $oldLines) líneas" -ForegroundColor $(if($newLines -gt $oldLines){'Green'}else{'Red'})
            Write-Host "Líneas eliminadas: $removed" -ForegroundColor Red
            Write-Host "Líneas añadidas: $added" -ForegroundColor Green
        }
        "4" {
            code -d "EP_Cesar\R\$file" "R\$file"
            Write-Host "Archivos abiertos en VS Code para comparación:" -ForegroundColor Green
            Write-Host "  IZQUIERDA (nuevo): EP_Cesar/R/$file" -ForegroundColor Cyan
            Write-Host "  DERECHA (actual): R/$file" -ForegroundColor Yellow
            Write-Host "  Puedes usar las flechas en el editor para mover cambios de izquierda a derecha" -ForegroundColor White
        }
        "s" {
            Write-Host "Saltando archivo...`n" -ForegroundColor Yellow
        }
        "q" {
            Write-Host "Saliendo...`n" -ForegroundColor Red
            return
        }
        default {
            Write-Host "Opción no válida`n" -ForegroundColor Red
            $index--
        }
    }
    
    Write-Host "`n========================================`n"
    $index++
    
    if ($index -le $changedFiles.Count) {
        $continue = Read-Host "Presiona Enter para continuar al siguiente archivo o 'q' para salir"
        if ($continue -eq 'q') {
            Write-Host "Saliendo...`n" -ForegroundColor Red
            return
        }
    }
}

Write-Host "¡Revisión completada!" -ForegroundColor Green
