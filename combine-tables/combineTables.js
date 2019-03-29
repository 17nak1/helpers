

fs= require('fs')
combineTables = function (modeltype, runs = [1]) {
  let allSet = [] 
  for ( run = 0; run < runs.length; run++){
    var table = [], dataset = []
    var temp
    var data = fs.readFileSync('./' + modeltype + '_' + runs[run] +'.csv').toString()
    var lines = data.split('\n')
    for (let i = 1; i < lines.length; i++) {
      temp = lines[i].split(',')
      if (temp[temp.length - 1] !== '"NaN"')
        table.push(temp)
    }
    table.sort(sortFunction)
    var hashMap = {}

    table.forEach(function(arr){
      // If your subArrays can be in any order, you can use .sort to have consistant order
      hashMap[arr.join("|")] = arr;
    });

    var dataset = Object.keys(hashMap).map(function(k){
      return hashMap[k]
    })
    console.log(run, dataset.length)
    allSet.push(...dataset) 
  }
  allSet.sort(sortFunction)
  const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['S', 'E', 'I', 'R', 'H'],
    path: './all.csv'
  })
   
  csvWriter.writeRecords(allSet)
    .then(() => {
    console.log('...done')
  })
}
combineTables ('DeterministicSEIR', ['mu', 'psi','rho'])
function sortFunction(a, b) {
  if (Number(a[a.length - 1]) === Number(b[a.length - 1])) {
    return 0
  }
  else {
    return (Number(a[a.length - 1]) < Number(b[a.length - 1])) ? 1 : -1;
  }
}