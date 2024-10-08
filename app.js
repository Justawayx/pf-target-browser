// NPM packages
const express = require('express')
const handlebars = require('express-handlebars')
const sqlite3 = require('sqlite3')
const fs = require('fs')
const https = require('https');
const http = require('http');
var Promise = require('promise');

// Load preprocessed JSON files
const global_RMSDs = JSON.parse(fs.readFileSync('data/global_RMSDs.json', 'utf8'))
const BDB_target_info = JSON.parse(fs.readFileSync('data/BDB_target_info.json', 'utf8'))
const Hs_TMalign_scores = JSON.parse(fs.readFileSync('data/Hs_TMalign_scores.json', 'utf8'))
const Hs_TMalign_seqIds = JSON.parse(fs.readFileSync('data/Hs_TMalign_seqIds.json', 'utf8'))
const Hs_TMalign_RMSDs = JSON.parse(fs.readFileSync('data/Hs_TMalign_RMSDs.json', 'utf8'))
const MolecularWeights = JSON.parse(fs.readFileSync('data/MolecularWeights.json', 'utf8'))
const ProteinLengths = JSON.parse(fs.readFileSync('data/ProteinLengths.json', 'utf8'))
const IsoelectricPoints = JSON.parse(fs.readFileSync('data/IsoelectricPoints.json', 'utf8'))

const ResistomeNumSampleswithInterestingMissenseMutations_list = JSON.parse(fs.readFileSync('data/ResistomeNumSampleswithInterestingMissenseMutations.json', 'utf8'))
const ResistomeNumSampleswithMissenseMutations_list = JSON.parse(fs.readFileSync('data/ResistomeNumSampleswithMissenseMutations.json', 'utf8'))
const ResistomeNumSampleswithDisruptiveMutations_list = JSON.parse(fs.readFileSync('data/ResistomeNumSampleswithDisruptiveMutations.json', 'utf8'))

const PlasmoDBTotalSNPs_list = JSON.parse(fs.readFileSync('data/PlasmoDBTotalSNPs.json', 'utf8'))
const PlasmoDBNoncodingSNPs_list = JSON.parse(fs.readFileSync('data/PlasmoDBNoncodingSNPs.json', 'utf8'))
const PlasmoDBSynonymousSNPs_list = JSON.parse(fs.readFileSync('data/PlasmoDBSynonymousSNPs.json', 'utf8'))
const PlasmoDBNonsynonymousSNPs_list = JSON.parse(fs.readFileSync('data/PlasmoDBNonsynonymousSNPs.json', 'utf8'))
const PlasmoDBStopCodonSNPs_list = JSON.parse(fs.readFileSync('data/PlasmoDBStopCodonSNPs.json', 'utf8'))

// SSL
const privateKey  = fs.readFileSync('/root/certs/selfsigned.key', 'utf8');
const certificate = fs.readFileSync('/root/certs/selfsigned.crt', 'utf8');
const credentials = {key: privateKey, cert: certificate};

// Start app
const app = express()
const httpServer = http.createServer(app);
const httpsServer = https.createServer(credentials, app);
const http_port = 80 // 3000
const https_port = 443 // 3000

// Set up Handlebars
app.set('view engine', 'handlebars')
app.engine('handlebars', handlebars.engine({
	layoutsDir: __dirname + '/views/layouts',
	helpers: {
		check: function(value) { return (!value) ? 'N/A' : value; },
		abbrcheck: function(value) { return (!value || value == "N/A") ? '' : ` (${value})`; },
		divide1000: function(value) { return value/1000; },
		round: function(value) { return Number(value).toFixed(2) },
		roundnodata: function(value) { return (!value) ? '<span style="color:gray;">nd</span>' : Number(value).toFixed(2) },
		textnodata: function(value) { return (!value) ? '<span style="color:gray;">nd</span>' : value },
		expr: function(value) { return (value == 0) ? `<em>${value}</em>` : `<em style="color:lime">${Number(value).toFixed(2)}</em>` },
	}
}))
app.use(express.static('public'))

// Load database
const db = new sqlite3.Database('pf_gene_anno.db')

function getPubData(geneID) {
	return new Promise((resolve, reject) => {
		db.all(`SELECT * FROM pmid WHERE GeneID = '${geneID}'`, (err, rows) => {
			if (err) reject(err);
			else resolve(rows);
		})})
}

function getSTRINGLink(req_id) {
	return new Promise((resolve, reject) => {
		https.get(`https://string-db.org/api/json/get_link?identifiers=${req_id}`, (resp) => {
			let data = '';
			resp.on('data', (chunk) => { data += chunk; });
			resp.on('end', () => { resolve(JSON.parse(data)[0]); });
		}).on("error", (err) => { reject(err); });
	})
}

async function getAnnoData(req_id) {
	const STRINGlink = await getSTRINGLink(req_id);
	return new Promise((resolve, reject) => {
		db.all(`SELECT * FROM anno WHERE GeneID = '${req_id}' OR TranscriptID = '${req_id}'`, (err, rows) => {
			if (err || rows.length == 0) reject(err);
			else resolve({"anno_rows": rows, "STRING_link": STRINGlink});
		})})
}

function getECData(req_id) {
	return new Promise((resolve, reject) => {
		db.all(`SELECT * FROM ec_inhibitors WHERE GeneID = '${req_id}'`, (err, rows) => {
			if (err) reject(err);
			else resolve(rows);
		})})
}

// Search bar

function search() {
	let input = document.getElementById('searchbar').value
	console.log(input)
}

// Helper functions
function zip(a, b) {
	c = a.map(function(a_elem, i) { return [a_elem, b[i]] })
	return c
}

function percentRank(array, n) {
	let i = 0.0
	for (let j = 0; j < array.length; j++) {
		if (n <= array[j]) break;
		i += 1.0
	}
	return (i/array.length)*100
}

function dynamicTableHeight(num_rows, row_height=22, extra=38) {
	if (num_rows > 4) return '';
	else return ` style="height:${(num_rows*row_height)+extra}px"`
}

function make_GO_table(cIDs, cs, fIDs, fs, pIDs, ps) {
	var GO_url_prefix = 'https://amigo.geneontology.org/amigo/term'
	var table_data = []
	
	table_data = table_data.concat(cIDs.map(function(x, i) { if (x != 'N/A') return ["Component", `<a href="${GO_url_prefix}/${x}" target="_blank">${x}</a>`, cs[i]] }).filter(function(x) { return x }))	
	table_data = table_data.concat(fIDs.map(function(x, i) { if (x != 'N/A') return ["Function", `<a href="${GO_url_prefix}/${x}" target="_blank">${x}</a>`, fs[i]] }).filter(function(x) { return x }))
	table_data = table_data.concat(pIDs.map(function(x, i) { if (x != 'N/A') return ["Process", `<a href="${GO_url_prefix}/${x}" target="_blank">${x}</a>`, ps[i]] }).filter(function(x) { return x }))
	if (table_data.length == 0) return null;
	else return table_data
}

function make_EC_inhibitors_table(ec_rows, num_inhib=12) {
	table_data = []
	let visual_row_count = 0
	
	if (ec_rows.length == 0) return [null, 0];
	for (let i=0; i < ec_rows.length; i++) {
		
		ec_number = ec_rows[i].ECNumber
		ec_elem = `<a href="https://www.brenda-enzymes.org/enzyme.php?ecno=${ec_number}" target=_blank>${ec_number}</a>`
		
		inhib_elem = ''
		if (!ec_rows[i].Inhibitors) {
			inhib_elem = "<i>No BRENDA inhibitors</i>";
			visual_row_count += 1;
		}
		else if (ec_rows[i].Inhibitors == 'BRENDA SOAP query unsuccessful') {
			inhib_elem = "<i>BRENDA SOAP query unsuccessful</i>";
			visual_row_count += 1;
		}
		else {
			inhibitors = ec_rows[i].Inhibitors.split(';')
			inhibitors_subset = []
			for (let j = 0; j < inhibitors.length; j++) {
				if (!((inhibitors[j].includes('+') || inhibitors[j].includes('-')) && (inhibitors[j].length <= 3) || (inhibitors[j].length == 4 && inhibitors[j][1].toLowerCase() == inhibitors[j][1]))) {
					inhibitors_subset.push(inhibitors[j])
				}
			}			
			if (inhibitors_subset.length > num_inhib) {
				for (let j=0; j < num_inhib; j++) {
					inhib_elem += `<span class="blob">${inhibitors_subset[j]}</span>`
					visual_row_count += 2;
				}
				inhib_elem += '<span class="blob">...</span>';
			} else {
				for (let j=0; j < inhibitors_subset.length; j++) {
					inhib_elem += `<span class="blob">${inhibitors_subset[j]}</span>`
				}
				if (inhibitors_subset.length < 6) visual_row_count += 1;
				else visual_row_count += 2;
			}
		}
		table_data.push([ec_elem, ec_rows[i].ECName, inhib_elem])		
	}
	return [table_data, visual_row_count]
}

function make_Pf7_table(JSON_str) {	
	const effect_list = ['synonymous', 'disruptive', 'missense']
	const freq_bin_list = ['common', 'rare', 'doubleton', 'singleton']
	table_data = []	
	if (JSON_str) {
		effect_freq_JSON = JSON.parse(JSON_str)
		for (let i=0; i < effect_list.length; i++) {
			effect = effect_list[i]
			row_name = effect
			row = freq_bin_list.map(function(freq_bin) { return (freq_bin in effect_freq_JSON[effect]) ? effect_freq_JSON[effect][freq_bin] : 0 })
			table_data.push([row_name].concat(row))
		}
	}
	return ((table_data.length == 0) ? null : table_data)
}

function make_BDB_table(bdbHOGENOM, bdbOMA, bdbOrthoDB, bdbOrthoMCL, bdbOrthoMCLBLAST) {
	var uniprot_source_dict = {}
	var ostrs = [bdbHOGENOM, bdbOMA, bdbOrthoDB, bdbOrthoMCL, bdbOrthoMCLBLAST]
	const source_list = ["HOGENOM", "OMA", "OrthoDB", "OrthoMCL", "OrthoMCL BLAST"]
	
	for (let i=0; i < ostrs.length; i++) {
		if (!ostrs[i]) continue;
		olist = ostrs[i].split(',')
		source = source_list[i]
		for (let j=0; j < olist.length; j++) {
			if (olist[j] == '') continue;
			uniprot = olist[j]
			if (!uniprot_source_dict[uniprot]) uniprot_source_dict[uniprot] = [];
			uniprot_source_dict[uniprot].push(source)
		}
	}
	
	table_data = []
	var all_uniprots = Object.keys(uniprot_source_dict)
	for (let i=0; i < all_uniprots.length; i++) {
		uniprot = all_uniprots[i]
		tinfo = BDB_target_info[uniprot]
		ligand_divs = tinfo.BDBLigandNames.map(function(name) { return `<span class="blob">${name}</span>` })
		if (ligand_divs.length > 4) ligand_divs = ligand_divs.slice(0, 4).concat(['<span class="blob">...</span>'])
		
		table_data.push([`<a href="${tinfo.BDBTargetLink}" target="_blank">${uniprot}</a>`, tinfo.BDBTargetName, tinfo.BDBTargetOrganism, `<span style="font-size:11px">${uniprot_source_dict[uniprot].join(', ')}</span>`, ligand_divs.join('')])
	}
	
	return table_data
}

function make_domain_table(iIDs, is, pIDs, ps, sIDs, ss) {
	var table_data = []
	table_data = table_data.concat(iIDs.map(function(x, i) { return ["InterPro", x, is[i]] }))
	table_data = table_data.concat(pIDs.map(function(x, i) { return ["PFam", x, ps[i]] }))
	table_data = table_data.concat(sIDs.map(function(x, i) { return ["Superfamily", x, ss[i]] }))
	return table_data
}

function make_genome_browser_link(loc_str) {
	var items = loc_str.split(':')
	var chrom = items[0]
	var suffix_items = items[1].split('(')[0].split('..')
	var url = `https://plasmodb.org/plasmo/app/jbrowse?data=%2Fa%2Fservice%2Fjbrowse%2Ftracks%2Fdefault&loc=${chrom}%3A${suffix_items[0]}..${suffix_items[1]}&tracks=gene%2CRNA-Seq%20Evidence%20for%20Introns%2CPfalciparum3D7%20combined%20RNAseq%20plot&highlight=`
	return `<a href="${url}" target="_blank">${loc_str}</a>`
}

function make_uniprots_str(uniprots) {
	var uniprot_links = uniprots.map(function(x) { return `<a href="https://www.uniprot.org/uniprotkb/${x}/entry" target="_blank">${x}</a>` })
	return uniprot_links.join(', ')
}

function make_pdbs_str(pdbs) {
	var pdb_links = pdbs.map(function(x) { return `<a href="https://www.rcsb.org/structure/${x}" target="_blank">${x}</a>` })
	return pdb_links.join(', ')
}

function getGreenToRed(percent){
	r = percent<50 ? 255 : Math.floor(255-(percent*2-100)*255/100);
	g = percent>50 ? 255 : Math.floor((percent*2)*255/100);
	return 'rgb('+r+','+g+',0)';
}

function make_alphafill_str(transplant_ID, transplant_name, local_RMSD, hit_PDB, global_RMSD) {
	
	var local_RMSD_color = 'white'
	local_RMSD = Number(local_RMSD).toFixed(2)
	if (local_RMSD > 4) local_RMSD_color = 'red';
	else if (local_RMSD > 1) local_RMSD_color = 'yellow';
	
	if (!transplant_name) transplant_name = ''; // No transplant full name
	
	global_RMSD = Number(global_RMSD).toFixed(2)
	var global_RMSD_color = getGreenToRed(100-percentRank(global_RMSDs, global_RMSD))
	
	if (!hit_PDB) return '<i>No AlphaFill hits</i>'
	else return `<a href="https://www.rcsb.org/ligand/${transplant_ID}" target=_blank>${transplant_ID}</a> (${transplant_name.split(';')[0].toLowerCase()}, Local RMSD=<em style="color:${local_RMSD_color}">${local_RMSD}</em>) with <a href="https://www.rcsb.org/structure/${hit_PDB}" target=_blank>${hit_PDB}</a> (Global RMSD=<em style="color:${global_RMSD_color}">${global_RMSD}</em>)`
}

function format_missense(muts, interesting_muts) {
	if (!muts) return "None";
	var mut_list = muts.split(',')
	var int_mut_list = (!interesting_muts) ? [] : interesting_muts.split(',')
	var formatted_muts = []
	for (let i = 0; i < mut_list.length; i++) {
		mut = mut_list[i]
		if (int_mut_list.includes(mut)) formatted_muts.push(`<b style="color:green">${mut}</b>`);
		else formatted_muts.push(mut);
	}
	return formatted_muts.join(', ')
}

function format_compounds(cmpds, interesting_cmpds) {
	if (!cmpds) return "None";
	var cmpd_list = cmpds.split(',')
	var int_cmpd_list = (!interesting_cmpds) ? [] : interesting_cmpds.split(',')
	var formatted_cmpds = []
	for (let i = 0; i < cmpd_list.length; i++) {
		cmpd = cmpd_list[i]
		if (int_cmpd_list.includes(cmpd)) formatted_cmpds.push(`<b style="color:green">${cmpd}</b>`);
		else formatted_cmpds.push(cmpd);
	}
	return formatted_cmpds.join(', ')
}

function format_EC(ec_str) {
	if (ec_str == 'N/A') return '<i>None</i>'
	var items = ec_str.split(';')
	for (let i = 0; i < items.length; i++) {
		ec_number = items[i].split(' (')[0]
		if (ec_number.includes('involved in cellular and subcellular movement.')) continue;
		ec_str = ec_str.replace(ec_number, `<a href="https://www.brenda-enzymes.org/enzyme.php?ecno=${ec_number}" target=_blank>${ec_number}</a>`)
	}
	return ec_str
}

var PlasmoGEM_phenotype_color_dict = {
	"Essential": 'green',
	"Slow": 'olive',
	'Dispensable': 'red',
	"Fast": 'pink',
	"Insufficient data": 'gray',
}

function get_MFS_color(MFS) {
	if (MFS < 0) return "lime";
	else return "red";
}

function get_RGR_color(relative_growth_rate) {	
	var prefix = relative_growth_rate.split(' ')[0]
	if (prefix == 'nan') return "white";
	else return getGreenToRed(100*(1-parseFloat(prefix)));	
}

function get_NumIns_color(numIns) {
	if (numIns == 0) return "lime";
	else if (numIns == 1) return "yellow";
	else return "red";
}

var Zhang_phenotype_color_dict = {
	"Mutable in CDS": 'red',
	"Non - Mutable in CDS": 'green',
}

var RMgmDB_phenotype_translation = {
	"nt": "Not tested",
	"nd": "Not different from wild type",
	"X": "Different from wild type",
	null: null,
}

var RMgmDB_phenotype_color_dict = {
	"nt": "gray",
	"nd": "red",
	"X": "green",
	null: null,
}

// Home
app.get('/', (req, res) => {
	res.redirect('/about');
})

// About
app.get('/about', (req, res) => {
	res.render('about', {layout : 'general_layout'});
})

// Gene list (serve static page for efficiency)
app.get('/genelist', (req, res) => {
	res.sendfile('public/genelist.html');
})

// Search

function check_search_match(query, data) {
	let lowercase_query = query.toLowerCase();
	for (let i= 0; i < data.length; i++) {
		if (data[i].toLowerCase().includes(lowercase_query)) return true;
	}
	return false;
}

app.get('/search', (req, res) => {
		
		var search_query = req.query.q
		
		db.all(`SELECT GeneID, GeneSymbol, ProductDescription FROM anno`, (err, rows) => {
		if (err) { 
			res.status(500).json({ error: err.message });
			return; 
		}
		var genes = [];
		var tups = [];
		for (let i = 0; i < rows.length; i++) {
			let search_match = check_search_match(search_query, [rows[i].GeneID, rows[i].GeneSymbol, rows[i].ProductDescription])
			if (search_match) {
				if (!genes.includes(rows[i].GeneID)) {
					genes.push(rows[i].GeneID);
					tups.push([`<a href="./${rows[i].GeneID}" target=_blank>${rows[i].GeneID}</a>`, (rows[i].GeneSymbol == 'N/A') ? '' : rows[i].GeneSymbol, rows[i].ProductDescription]);
				}
			}
		}
		if (genes.length == 1) {
			res.redirect(`/${genes[0]}`);
		}
		else {
			res.render('search', {layout: 'general_layout', data: {'geneData': tups, 'query': search_query}});
		}
	})
})

// Individual gene pages
app.get('/PF3D7*', (req, res) => {
	
	var new_layout = "genepage_layout"
	var new_main = "main2"
	
	console.log(req.url)
	var req_id = req.url.split('?')[0].slice(1,)
	
	getAnnoData(req_id).then(obj => {	
		return obj
	}).then(obj => {
		getECData(req_id).then(ec_rows => {
			return {'anno': obj['anno_rows'], 'STRING_link': obj['STRING_link'], 'ec': ec_rows}
		}).then(data => {
			
		getPubData(req_id).then(pmid_rows => {
			
			anno_rows = data['anno']; ec_rows = data['ec']; STRING_link = data['STRING_link']
			// anno_rows: from anno table
			// pmid_rows: from pmid table
			if (anno_rows.length == 0) {
				res.writeHead(400, { "Content-Type": "text/plain" });
				res.end("Invalid Gene/Transcript ID\n");
				return;
			}
			
			// Use anno table data as base
			var row = anno_rows[0]
			
			// Reformat pmid table data
			
			var pmidTableData = []
			for (let i = 0; i < pmid_rows.length; i++) {
				pmid_row = pmid_rows[i]
				doiElement = ((!pmid_row.DOI) || (!pmid_row.DOI.includes('/'))) ? '' : `<a href="https://doi.org/${pmid_row.DOI}" target=_blank>${pmid_row.DOI}</a>`
				authors = (!pmid_row.Authors) ? [] : pmid_row.Authors.split(',')
				authorElement = (authors.length > 4) ? [authors[0], authors[1], '...', authors[authors.length-1]].join(', ') : authors.join(', ')				
				pmidTableData.push([pmid_row.PMID, pmid_row.Title, authorElement, doiElement])
			}
			
			row.pmidTableData = (pmidTableData.length == 0) ? null : pmidTableData
			row.STRINGlink = STRING_link
			
			// Reformat anno table data
			row.ECNumbersFormatted = format_EC(row.ECNumbers)
			let [ecTableData, visual_row_count] = make_EC_inhibitors_table(ec_rows, (new_layout == 'index2') ? 4 : 12)
			
			row.ECTableData = ecTableData
			if (row.ECTableData) {
				row.ECTableHeight = dynamicTableHeight(visual_row_count)
			}
			row.HsGeneTable = (!row.HsGeneIDs) ? null : zip(row.HsGeneIDs.split(';'), row.HsGeneDescriptions.split(';'))
			if (row.HsGeneTable) {
				row.HsGeneTableHeight = dynamicTableHeight(row.HsGeneTable.length)
			}
			row.GOTableData = make_GO_table(row.CuratedGOComponentIDs.split(';'), row.CuratedGOComponents.split(';'), row.CuratedGOFunctionIDs.split(';'), row.CuratedGOFunctions.split(';'), row.CuratedGOProcessIDs.split(';'), row.CuratedGOProcesses.split(';'))
			if (row.GOTableData) {
				row.GOTableHeight = dynamicTableHeight(row.GOTableData.length)
			}
			row.ProteinDomainData = make_domain_table(row.InterproID.split(';'), row.InterproDescription.split(';'), row.PFamID.split(';'), row.PFamDescription.split(';'), row.SuperfamilyID.split(';'), row.SuperfamilyDescription.split(';'))
			row.GenomicLocationElement = make_genome_browser_link(row.GenomicLocation)			
			row.UniprotIDLinks = (row.UniprotIDs == 'N/A') ? '<i>None</i>' : make_uniprots_str(row.UniprotIDs.split(','))
			row.PDBIDLinks = (row.PDBIDs) ? make_pdbs_str(row.PDBIDs.split(',')) : '<i>None</i>'
			row.AlphaFillSummary = make_alphafill_str(row.AlphaFillHitTransplantID, row.AlphaFillHitTransplantName, row.AlphaFillHitTransplantLocalRMSD, row.AlphaFillHitPDBID, row.AlphaFillHitGlobalRMSD)
			row.Pf7TableHomoData = make_Pf7_table(row.Pf7NumberVariantsbyEffectxHomozygousPrevalenceJSON)
			row.Pf7TableAnyData = make_Pf7_table(row.Pf7NumberVariantsbyEffectxAnyPrevalenceJSON)
			row.BindingDBTableData = make_BDB_table(row.BindingDBMatchHOGENOM, row.BindingDBMatchOMA, row.BindingDBMatchOrthoDB, row.BindingDBMatchOrthoMCL, row.BindingDBMatchOrthoMCLBLAST)
			if (row.BindingDBTableData) {
				row.BindingDBTableHeight = dynamicTableHeight(row.BindingDBTableData.length, 30)
			}
			row.RMgmDBPhenotype = RMgmDB_phenotype_translation[row.RMgmDBABSDifferentfromWT]
			row.RMgmDBPhenotypeColor = RMgmDB_phenotype_color_dict[row.RMgmDBABSDifferentfromWT]			
			row.FormattedMissenseMutations = format_missense(row.ResistomeMissenseMutations, row.ResistomeInterestingMissenseMutations)
			row.FormattedResistomeCompounds = format_compounds(row.ResistomeCompoundsofSampleswithMissenseMutations, row.ResistomeCompoundsofSampleswithInterestingMissenseMutations)			
			row.ResistomeNumSampleswithMissenseMutations = row['Resistome#SampleswithMissenseMutations']
			row.ResistomeNumSampleswithInterestingMissenseMutations = row['Resistome#SampleswithInterestingMissenseMutations']
			row.ResistomeNumSampleswithDisruptiveMutations = row['Resistome#SampleswithDisruptiveMutations']
			
			// Colors
			row.HsTMalignScoreColor = getGreenToRed(percentRank(Hs_TMalign_scores, row.HsTMalignScore))
			row.HsTMalignRMSDColor = getGreenToRed(100-percentRank(Hs_TMalign_RMSDs, row.HsTMalignRMSD))
			row.HsTMalignSeqIdentityColor = getGreenToRed(percentRank(Hs_TMalign_seqIds, row.HsTMalignSeqIdentity))
			
			row.PlasmoGEMPhenotypeColor = (row.PlasmoGEMPhenotype) ? PlasmoGEM_phenotype_color_dict[row.PlasmoGEMPhenotype] : "black"
			if (row.PlasmoGEMRelativeGrowthRate) {
				row.PlasmoGEMRelativeGrowthRateColor = get_RGR_color(row.PlasmoGEMRelativeGrowthRate)
			}
			
			row.ZhangPhenotypeColor = (row.ZhangPhenotype) ? Zhang_phenotype_color_dict[row.ZhangPhenotype] : "black"
			row.ZhangMISColor = getGreenToRed(100*(1-row.ZhangOriginalMIS))
			row.ZhangMFSColor = get_MFS_color(row.ZhangOriginalMFS)
			row.ZhangNumInsertionsinCDSColor = get_NumIns_color(row.ZhangNumInsertionsinCDS)
			
			row.ProteinLengthColor = getGreenToRed(percentRank(ProteinLengths, row.ProteinLength))
			row.MolecularWeightColor = getGreenToRed(percentRank(MolecularWeights, row.MolecularWeight))
			row.IsoelectricPointColor = getGreenToRed(percentRank(IsoelectricPoints, row.IsoelectricPoint))
			
			row.PlasmoDBTotalSNPsColor = getGreenToRed(100-percentRank(PlasmoDBTotalSNPs_list, row.PlasmoDBTotalSNPs))
			row.PlasmoDBNoncodingSNPsColor = getGreenToRed(100-percentRank(PlasmoDBNoncodingSNPs_list, row.PlasmoDBNoncodingSNPs))
			row.PlasmoDBSynonymousSNPsColor = getGreenToRed(100-percentRank(PlasmoDBSynonymousSNPs_list, row.PlasmoDBSynonymousSNPs))
			row.PlasmoDBNonsynonymousSNPsColor = getGreenToRed(100-percentRank(PlasmoDBNonsynonymousSNPs_list, row.PlasmoDBNonsynonymousSNPs))
			row.PlasmoDBStopCodonSNPsColor = getGreenToRed(100-percentRank(PlasmoDBStopCodonSNPs_list, row.PlasmoDBStopCodonSNPs))
			
			row.ResistomeNumSampleswithInterestingMissenseMutationsColor = getGreenToRed(percentRank(ResistomeNumSampleswithInterestingMissenseMutations_list, row.ResistomeNumSampleswithInterestingMissenseMutations))
			row.ResistomeNumSampleswithMissenseMutationsColor = getGreenToRed(percentRank(ResistomeNumSampleswithMissenseMutations_list, row.ResistomeNumSampleswithMissenseMutations))
			row.ResistomeNumSampleswithDisruptiveMutationsColor = getGreenToRed(percentRank(ResistomeNumSampleswithDisruptiveMutations_list, row.ResistomeNumSampleswithDisruptiveMutations))
			
			// Flags
			row.HasHumanOrtholog = !(!row.HsMostSimilarUniprotID)
			row.HasBindingDBOrtholog = (row.BindingDBMatchOrthoMCLBLAST || row.BindingDBMatchAnyExact)
			row.HasZhang = !(!row.ZhangOriginalMIS)
			
			res.render(new_main, {layout : new_layout, data : row});	
			
		}).catch(err => { 
			res.status(500).json({ error: err.message });
			return; 
		})
	}).catch(err => { 
			res.status(500).json({ error: err.message });
			return; 
		})
	}).catch(err => { 
			res.status(500).json({ error: err.message });
			return; 
		})
})

//httpServer.listen(http_port);
//httpsServer.listen(https_port);
httpServer.listen(80);
httpsServer.listen(443);
