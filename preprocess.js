const sqlite3 = require('sqlite3')
const fs = require('fs')

const db = new sqlite3.Database('pf_gene_anno.db')

db.all(`SELECT * FROM anno`, (err, rows) => {
		if (err) { 
			console.log( err.message );
			return; 
		}
		
		var global_RMSDs = []
		for (let i = 0; i < rows.length; i++) {
			if (rows[i].AlphaFillHitGlobalRMSD) {
				global_RMSDs.push(parseFloat(rows[i].AlphaFillHitGlobalRMSD));
			}
		}
		fs.writeFile("data/global_RMSDs.json", JSON.stringify(global_RMSDs.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		
		var Hs_TMalign_scores = []
		var Hs_TMalign_seqIds = []
		var Hs_TMalign_RMSDs = []
		
		for (let i = 0; i < rows.length; i++) {
			if (rows[i].HsTMalignScore) {
				Hs_TMalign_scores.push(parseFloat(rows[i].HsTMalignScore));
				Hs_TMalign_seqIds.push(parseFloat(rows[i].HsTMalignSeqIdentity));
				Hs_TMalign_RMSDs.push(parseFloat(rows[i].HsTMalignRMSD));
			}
		}
		fs.writeFile("data/Hs_TMalign_scores.json", JSON.stringify(Hs_TMalign_scores.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/Hs_TMalign_seqIds.json", JSON.stringify(Hs_TMalign_seqIds.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/Hs_TMalign_RMSDs.json", JSON.stringify(Hs_TMalign_RMSDs.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		
		var ProteinLengths = []
		var MolecularWeights = []
		var IsoelectricPoints = []
		
		for (let i = 0; i < rows.length; i++) {
			ProteinLengths.push(parseInt(rows[i].ProteinLength));
			MolecularWeights.push(parseFloat(rows[i].MolecularWeight));
			IsoelectricPoints.push(parseFloat(rows[i].IsoelectricPoint));
		}
		fs.writeFile("data/ProteinLengths.json", JSON.stringify(ProteinLengths.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/MolecularWeights.json", JSON.stringify(MolecularWeights.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/IsoelectricPoints.json", JSON.stringify(IsoelectricPoints.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
			
		var ResistomeNumSampleswithInterestingMissenseMutations = []
		var ResistomeNumSampleswithMissenseMutations = []
		var ResistomeNumSampleswithDisruptiveMutations = []
		
		for (let i = 0; i < rows.length; i++) {
			if (rows[i]['Resistome#SampleswithMissenseMutations']) {
				ResistomeNumSampleswithInterestingMissenseMutations.push(parseInt(rows[i]['Resistome#SampleswithInterestingMissenseMutations']));
				ResistomeNumSampleswithMissenseMutations.push(parseInt(rows[i]['Resistome#SampleswithMissenseMutations']));
				ResistomeNumSampleswithDisruptiveMutations.push(parseInt(rows[i]['Resistome#SampleswithDisruptiveMutations']));
			}
		}
		fs.writeFile("data/ResistomeNumSampleswithInterestingMissenseMutations.json", JSON.stringify(ResistomeNumSampleswithInterestingMissenseMutations.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/ResistomeNumSampleswithMissenseMutations.json", JSON.stringify(ResistomeNumSampleswithMissenseMutations.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/ResistomeNumSampleswithDisruptiveMutations.json", JSON.stringify(ResistomeNumSampleswithDisruptiveMutations.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		
		var PlasmoDBTotalSNPs = []
		var PlasmoDBNoncodingSNPs = []
		var PlasmoDBSynonymousSNPs = []
		var PlasmoDBNonsynonymousSNPs = []
		var PlasmoDBStopCodonSNPs = []
		
		for (let i = 0; i < rows.length; i++) {
			PlasmoDBTotalSNPs.push(parseInt(rows[i].PlasmoDBTotalSNPs));
			PlasmoDBNoncodingSNPs.push(parseInt(rows[i].PlasmoDBNoncodingSNPs));
			PlasmoDBSynonymousSNPs.push(parseInt(rows[i].PlasmoDBSynonymousSNPs));
			PlasmoDBNonsynonymousSNPs.push(parseInt(rows[i].PlasmoDBNonsynonymousSNPs));
			PlasmoDBStopCodonSNPs.push(parseInt(rows[i].PlasmoDBStopCodonSNPs));
		}
		fs.writeFile("data/PlasmoDBTotalSNPs.json", JSON.stringify(PlasmoDBTotalSNPs.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/PlasmoDBNoncodingSNPs.json", JSON.stringify(PlasmoDBNoncodingSNPs.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/PlasmoDBSynonymousSNPs.json", JSON.stringify(PlasmoDBSynonymousSNPs.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/PlasmoDBNonsynonymousSNPs.json", JSON.stringify(PlasmoDBNonsynonymousSNPs.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
		fs.writeFile("data/PlasmoDBStopCodonSNPs.json", JSON.stringify(PlasmoDBStopCodonSNPs.sort(function(a,b) { return a-b; })), function(err) { if (err) { console.log(err); } })
})



db.all(`SELECT * FROM bdb_targets`, (err, rows) => {
		if (err) { 
			console.log( err.message );
			return; 
		}
		var BDB_target_info = {}
		for (let i = 0; i < rows.length; i++) {
			var row = rows[i]
			BDB_target_info[row.Uniprot] = {
				BDBTargetName: row.BDBTargetName,
				BDBTargetOrganism: row.BDBTargetOrganism,
				BDBLigandNames: row.BDBLigandNames.split('; '),
				BDBTargetLink: row.BDBTargetLink
			}
		}
		fs.writeFile("data/BDB_target_info.json", JSON.stringify(BDB_target_info), function(err) { if (err) { console.log(err); } })
})