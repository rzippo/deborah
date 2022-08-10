$scenarios = Get-ChildItem -Recurse -Path "../scenarios/" -Filter "*.conf";

$commands = @(
    "lb",
    "lb --lb-experimental",
    "ludb",
    "ludb --ludb-nnested-sta",
    "foi-output-curve"
);

foreach($command in $commands)
{
    New-Item -Force -ItemType Directory $command;

    foreach($scenario in $scenarios)
    {
        $logName = "$($scenario.BaseName).log";
        $expression = "deborah $scenario --det-output --$command > '$command/$logName'";
        Write-Host $expression;
        Invoke-Expression $expression;
    }    
}