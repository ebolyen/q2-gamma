import loadWASM from './test.rs';

loadWASM().then(module => {
    const add = module.instance.exports['add'];
    console.log('foo:', add(1, 2));
})
